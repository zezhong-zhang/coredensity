"""
Phonon-Averaged Diffraction Analysis
====================================
Simple class-oriented approach for phonon configuration averaging.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import time
from typing import List, Dict, Tuple, Optional

import abtem
from abtem import Potential
from abtem.atoms import orthogonalize_cell
from scipy.interpolate import RegularGridInterpolator


class PhononDiffraction:
    """
    Class for phonon-averaged diffraction analysis for IAM and IAM+bonding.
    
    Workflow:
    1. Load configurations
    2. Sum potentials (IAM + bonding)  
    3. Run multislice simulations
    4. Visualize results
    """
    
    def __init__(self, data_dir: str = "pot-diff/pot_diff/data", 
                 energy: float = 100e3, max_angle: float = 50):
        """
        Initialize phonon diffraction analyzer.
        
        Parameters
        ----------
        data_dir : str
            Directory containing configuration files
        energy : float
            Electron beam energy in eV
        max_angle : float
            Maximum diffraction angle in mrad
        """
        self.data_dir = Path(data_dir)
        self.energy = energy
        self.max_angle = max_angle
        
        # Storage for results
        self.config_files: List[Path] = []
        
        # Store complex waves for proper elastic/inelastic separation
        self.iam_waves: List[np.ndarray] = []  # Complex waves
        self.total_waves: List[np.ndarray] = []  # Complex waves
        self.iam_intensities: List[np.ndarray] = []  # |wave|²
        self.total_intensities: List[np.ndarray] = []  # |wave|²
        
        # Elastic scattering (wave average)
        self.iam_elastic: Optional[np.ndarray] = None
        self.dft_elastic: Optional[np.ndarray] = None
        
        # Total scattering (intensity average)
        self.iam_total_intensity: Optional[np.ndarray] = None
        self.dft_total_intensity: Optional[np.ndarray] = None
        
        # Phonon scattering (Total - Elastic)
        self.iam_phonon: Optional[np.ndarray] = None
        self.total_phonon: Optional[np.ndarray] = None
        
        # abTEM objects for visualization and proper calculations
        self.iam_diffraction_template: Optional[abtem.DiffractionPatterns] = None
        self.iam_wave_template: Optional = None
        self.atoms_ortho: Optional = None
        
        # Metadata
        self.sampling: Optional[np.ndarray] = None
        self.extent: Optional[List[float]] = None
    
    def load_configurations(self, max_configs: Optional[int] = None) -> List[Path]:
        """Load configuration file paths."""
        pattern = str(self.data_dir / "*.pkl")
        self.config_files = sorted([Path(f) for f in glob.glob(pattern)])
        
        if max_configs:
            self.config_files = self.config_files[:max_configs]
        
        print(f"Found {len(self.config_files)} configurations")
        return self.config_files
    
    def _load_and_process_config(self, filepath: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Load single configuration and return complex waves and intensities.
        
        Returns
        -------
        tuple
            (iam_wave, total_wave, iam_intensity, total_intensity)
        """
        # Load DFT data
        data, atoms, _ = np.load(filepath, allow_pickle=True)
        
        # Z-integrate
        cell = atoms.get_cell()
        nz = data.shape[2]
        dz = np.linalg.norm(cell[2]) / nz
        dft_2d = np.sum(data, axis=2) * dz
        
        # Transform to orthogonal
        atoms_ortho = orthogonalize_cell(atoms, tolerance=0)
        dft_ortho, sampling = self._transform_to_orthogonal(dft_2d, atoms, atoms_ortho)
        
        # Create IAM potential
        iam_potential = Potential(
            atoms_ortho,
            parametrization="lobato",
            sampling=sampling,
            gpts=dft_ortho.shape,
            slice_thickness=10,
            projection="finite"
        ).build()
        
        # Create total potential (IAM + DFT bonding)
        total_potential = iam_potential.copy()
        total_potential.array[:] = np.array(iam_potential.array) + dft_ortho
        
        # Calculate exit waves (complex) and store abTEM objects
        planewave = abtem.PlaneWave(energy=self.energy)
        
        # Get complex exit waves as abTEM objects (not just arrays)
        iam_wave_obj = planewave.multislice(iam_potential).compute()
        total_wave_obj = planewave.multislice(total_potential).compute()
        
        # Calculate diffraction patterns from waves
        iam_diffraction = iam_wave_obj.diffraction_patterns(block_direct=True, max_angle=self.max_angle).compute()
        total_diffraction = total_wave_obj.diffraction_patterns(block_direct=True, max_angle=self.max_angle).compute()
        
        # Store metadata and templates from first successful config
        if self.sampling is None:
            self.sampling = np.array(iam_diffraction.sampling)
            # Convert extent to matplotlib format [xmin, xmax, ymin, ymax]
            extent = iam_diffraction.extent
            self.extent = [-extent[0]/2, extent[0]/2, -extent[1]/2, extent[1]/2]
            # Store templates for proper abTEM calculations
            self.iam_diffraction_template = iam_diffraction
            self.iam_wave_template = iam_wave_obj
            self.atoms_ortho = atoms_ortho
        
        # Extract arrays for storage
        iam_wave_array = np.array(iam_wave_obj.array)
        total_wave_array = np.array(total_wave_obj.array)
        iam_intensity = np.array(iam_diffraction.array)
        total_intensity = np.array(total_diffraction.array)
        
        return iam_wave_array, total_wave_array, iam_intensity, total_intensity
    
    def _transform_to_orthogonal(self, data_2d: np.ndarray, atoms, atoms_ortho) -> Tuple[np.ndarray, np.ndarray]:
        """Transform oblique data to orthogonal grid."""
        cell = atoms.get_cell()
        cell_ortho = atoms_ortho.get_cell()
        nx, ny = data_2d.shape
        
        # Calculate uniform sampling
        sampling_a = np.linalg.norm(cell[0, :2]) / nx
        sampling_b = np.linalg.norm(cell[1, :2]) / ny
        sampling_uniform = min(sampling_a, sampling_b)
        
        # Orthogonal grid dimensions
        nx_ortho = int(np.ceil(cell_ortho[0, 0] / sampling_uniform))
        ny_ortho = int(np.ceil(cell_ortho[1, 1] / sampling_uniform))
        
        # Tile and interpolate (simplified version)
        tile_size = 3
        data_tiled = np.tile(data_2d, (tile_size, tile_size))
        
        # Create interpolator
        frac_coords = (
            np.linspace(0, tile_size, nx * tile_size, endpoint=False),
            np.linspace(0, tile_size, ny * tile_size, endpoint=False)
        )
        interpolator = RegularGridInterpolator(frac_coords, data_tiled, 
                                             method='linear', bounds_error=False, fill_value=0)
        
        # Create orthogonal grid and transform coordinates
        x_ortho = np.arange(nx_ortho) * sampling_uniform
        y_ortho = np.arange(ny_ortho) * sampling_uniform
        X_ortho, Y_ortho = np.meshgrid(x_ortho, y_ortho, indexing='ij')
        
        # Transform to fractional coordinates
        cell_matrix = np.array([[cell[0, 0], cell[1, 0]], [cell[0, 1], cell[1, 1]]])
        cell_inv = np.linalg.inv(cell_matrix)
        
        ortho_points = np.column_stack([X_ortho.ravel(), Y_ortho.ravel()])
        frac_points = ortho_points @ cell_inv.T + tile_size // 2  # Center in tiled data
        
        # Interpolate
        data_ortho = interpolator(frac_points).reshape(nx_ortho, ny_ortho)
        
        return data_ortho, np.array([sampling_uniform, sampling_uniform])
    
    def process_all_configurations(self, max_configs: Optional[int] = None) -> Dict:
        """
        Process all configurations and calculate diffraction patterns.
        
        Returns
        -------
        dict
            Processing summary
        """
        if not self.config_files:
            self.load_configurations(max_configs)
        
        print(f"Processing {len(self.config_files)} configurations...")
        
        successful = 0
        failed = 0
        
        start_time = time.time()
        
        for i, filepath in enumerate(self.config_files, 1):
            try:
                print(f"  {i}/{len(self.config_files)}: {filepath.stem}")
                
                iam_wave, total_wave, iam_intensity, total_intensity = self._load_and_process_config(filepath)
                
                print(f"    Wave shapes: {iam_wave.shape}, Intensity shapes: {iam_intensity.shape}")
                
                self.iam_waves.append(iam_wave)
                self.total_waves.append(total_wave)
                self.iam_intensities.append(iam_intensity)
                self.total_intensities.append(total_intensity)
                successful += 1
                
            except Exception as e:
                print(f"    Error: {e}")
                failed += 1
        
        processing_time = time.time() - start_time
        
        # Calculate averages
        if self.iam_waves:
            self._calculate_elastic_inelastic_averages()
        
        summary = {
            'total': len(self.config_files),
            'successful': successful,
            'failed': failed,
            'processing_time': processing_time
        }
        
        print(f"Processing complete: {successful} successful, {failed} failed")
        print(f"Time: {processing_time:.1f} seconds")
        
        return summary
    
    def _calculate_elastic_inelastic_averages(self):
        """Calculate elastic and inelastic scattering components."""
        print("Calculating elastic and inelastic scattering...")
        
        # Check and ensure consistent shapes
        self._ensure_consistent_shapes()
        
        # Stack waves and intensities
        iam_wave_stack = np.stack(self.iam_waves, axis=0)
        total_wave_stack = np.stack(self.total_waves, axis=0)
        iam_intensity_stack = np.stack(self.iam_intensities, axis=0)
        total_intensity_stack = np.stack(self.total_intensities, axis=0)
        
        # Elastic scattering: |<wave>|² (coherent average)
        iam_wave_avg = np.mean(iam_wave_stack, axis=0)
        dft_wave_avg = np.mean(total_wave_stack, axis=0)
        
        # Convert waves to diffraction patterns for elastic scattering using abTEM
        self.iam_elastic = self._calculate_elastic_diffraction(iam_wave_avg)
        self.dft_elastic = self._calculate_elastic_diffraction(dft_wave_avg)
        
        # Total scattering: <|wave|²> (incoherent average)
        self.iam_total_intensity = np.mean(iam_intensity_stack, axis=0)
        self.dft_total_intensity = np.mean(total_intensity_stack, axis=0)
        
        # Phonon scattering: Total - Elastic
        self.iam_phonon = self.iam_total_intensity - self.iam_elastic
        self.total_phonon = self.dft_total_intensity - self.dft_elastic
        
        print(f"Processed {len(self.iam_waves)} configurations")
        
        # Debug: Check intensity scales
        print(f"Intensity scales:")
        print(f"  IAM total intensity: {np.sum(self.iam_total_intensity):.2e}")
        print(f"  IAM elastic intensity: {np.sum(self.iam_elastic):.2e}")
        print(f"  DFT total intensity: {np.sum(self.dft_total_intensity):.2e}")
        print(f"  DFT elastic intensity: {np.sum(self.dft_elastic):.2e}")
        
        print(f"Elastic vs Total intensity ratios:")
        print(f"  IAM: {np.sum(self.iam_elastic)/np.sum(self.iam_total_intensity):.4f}")
        print(f"  DFT: {np.sum(self.dft_elastic)/np.sum(self.dft_total_intensity):.4f}")
        
        # Check for unreasonable values
        if np.sum(self.iam_elastic) > np.sum(self.iam_total_intensity) * 1.1:
            print("WARNING: IAM elastic intensity is higher than total - normalization issue!")
        if np.sum(self.dft_elastic) > np.sum(self.dft_total_intensity) * 1.1:
            print("WARNING: DFT elastic intensity is higher than total - normalization issue!")
    
    def _wave_to_diffraction_intensity(self, wave_array: np.ndarray) -> np.ndarray:
        """
        Convert exit wave to diffraction pattern intensity with proper normalization.
        
        Parameters
        ----------
        wave_array : np.ndarray
            Complex exit wave array
            
        Returns
        -------
        np.ndarray
            Diffraction pattern intensity (properly normalized)
        """
        # FFT to get diffraction pattern
        fft_wave = np.fft.fftshift(np.fft.fft2(wave_array))
        
        # Calculate intensity |FFT|²
        intensity = np.abs(fft_wave)**2
        
        # Normalize by array size to match abTEM convention
        # abTEM normalizes FFT by the number of pixels
        intensity = intensity / (wave_array.size)
        
        # Ensure intensity matches the expected diffraction pattern size
        if hasattr(self, 'iam_intensities') and len(self.iam_intensities) > 0:
            target_shape = self.iam_intensities[0].shape
            
            if intensity.shape != target_shape:
                # Resize using center crop or pad
                def resize_array(arr, target_shape):
                    current_shape = arr.shape
                    
                    # Calculate crop/pad for each dimension
                    slices = []
                    for i, (current, target) in enumerate(zip(current_shape, target_shape)):
                        if current > target:
                            # Crop from center
                            start = (current - target) // 2
                            slices.append(slice(start, start + target))
                        elif current < target:
                            # This shouldn't happen with FFT, but handle it
                            slices.append(slice(None))
                        else:
                            slices.append(slice(None))
                    
                    cropped = arr[tuple(slices)]
                    
                    # If still different size, use simple resize
                    if cropped.shape != target_shape:
                        from scipy.ndimage import zoom
                        zoom_factors = [t/c for c, t in zip(cropped.shape, target_shape)]
                        cropped = zoom(cropped, zoom_factors, order=1)
                    
                    return cropped
                
                intensity = resize_array(intensity, target_shape)
        
        # Additional normalization to match the scale of original diffraction patterns
        if hasattr(self, 'iam_intensities') and len(self.iam_intensities) > 0:
            # Use the first intensity pattern as reference for scaling
            reference_intensity = self.iam_intensities[0]
            
            # Scale to match the total intensity of the reference
            current_total = np.sum(intensity)
            reference_total = np.sum(reference_intensity)
            
            if current_total > 0:
                intensity = intensity * (reference_total / current_total)
        
        return intensity
    
    def _calculate_elastic_diffraction(self, wave_avg: np.ndarray) -> np.ndarray:
        """
        Calculate elastic diffraction pattern from averaged wave using abTEM's methods.
        
        Parameters
        ----------
        wave_avg : np.ndarray
            Averaged complex exit wave
            
        Returns
        -------
        np.ndarray
            Properly normalized diffraction pattern intensity
        """
        if hasattr(self, 'iam_wave_template') and self.iam_wave_template is not None:
            # Create a new abTEM wave object with the averaged wave data
            elastic_wave = self.iam_wave_template.copy()
            
            # Handle shape mismatch by cropping or padding as needed
            if elastic_wave.array.shape != wave_avg.shape:
                print(f"Shape mismatch detected: template {elastic_wave.array.shape}, wave_avg {wave_avg.shape}")
                
                def resize_to_match(array, target_shape):
                    """Resize array to target shape by cropping or padding."""
                    current_shape = array.shape
                    
                    # If shapes are the same, return as is
                    if current_shape == target_shape:
                        return array
                    
                    # Handle each dimension
                    result = array
                    for dim in range(len(current_shape)):
                        current_size = result.shape[dim]
                        target_size = target_shape[dim]
                        
                        if current_size > target_size:
                            # Crop from center
                            start = (current_size - target_size) // 2
                            slices = [slice(None)] * len(result.shape)
                            slices[dim] = slice(start, start + target_size)
                            result = result[tuple(slices)]
                        elif current_size < target_size:
                            # Pad with zeros
                            pad_width = [(0, 0)] * len(result.shape)
                            pad_before = (target_size - current_size) // 2
                            pad_after = target_size - current_size - pad_before
                            pad_width[dim] = (pad_before, pad_after)
                            result = np.pad(result, pad_width, mode='constant', constant_values=0)
                    
                    return result
                
                # If wave_avg is smaller, pad it to match template
                if any(w < t for w, t in zip(wave_avg.shape, elastic_wave.array.shape)):
                    print(f"Padding wave_avg from {wave_avg.shape} to {elastic_wave.array.shape}")
                    wave_avg = resize_to_match(wave_avg, elastic_wave.array.shape)
                else:
                    # Otherwise resize the template to match wave_avg
                    print(f"Resizing template from {elastic_wave.array.shape} to {wave_avg.shape}")
                    elastic_wave.array = resize_to_match(elastic_wave.array, wave_avg.shape)
            
            elastic_wave.array[:] = wave_avg
            
            # Use abTEM's built-in diffraction calculation for proper normalization
            elastic_diffraction = elastic_wave.diffraction_patterns(
                block_direct=True, 
                max_angle=self.max_angle
            ).compute()
            
            return np.array(elastic_diffraction.array)
        
        else:
            # Fallback: manual calculation with better normalization
            print("Warning: Using fallback elastic calculation")
            
            # Calculate FFT (this is what abTEM does internally)
            fft_wave = np.fft.fft2(wave_avg)
            fft_wave = np.fft.fftshift(fft_wave)
            
            # Calculate intensity
            intensity = np.abs(fft_wave)**2
            
            # Normalize properly - this is the key fix!
            # abTEM normalizes by the sampling area, not the array size
            if hasattr(self, 'sampling') and self.sampling is not None:
                # Sampling is in Å/pixel, so sampling^2 is the pixel area in Å²
                pixel_area = self.sampling[0] * self.sampling[1]
                intensity = intensity * pixel_area
            
            # Match the scale of the original diffraction patterns
            if hasattr(self, 'iam_intensities') and len(self.iam_intensities) > 0:
                # Crop to match target shape first
                target_shape = self.iam_intensities[0].shape
                
                if intensity.shape != target_shape:
                    # Center crop
                    h, w = intensity.shape
                    th, tw = target_shape
                    
                    if h >= th and w >= tw:
                        start_h = (h - th) // 2
                        start_w = (w - tw) // 2
                        intensity = intensity[start_h:start_h+th, start_w:start_w+tw]
            
            return intensity
    
    def _ensure_consistent_shapes(self):
        """Ensure all arrays have consistent shapes by cropping to minimum size."""
        if not self.iam_waves:
            return
        
        # Find minimum dimensions across all arrays
        min_shape_waves = None
        min_shape_intensities = None
        
        # Check wave shapes
        for wave in self.iam_waves + self.total_waves:
            if min_shape_waves is None:
                min_shape_waves = wave.shape
            else:
                min_shape_waves = tuple(min(a, b) for a, b in zip(min_shape_waves, wave.shape))
        
        # Check intensity shapes  
        for intensity in self.iam_intensities + self.total_intensities:
            if min_shape_intensities is None:
                min_shape_intensities = intensity.shape
            else:
                min_shape_intensities = tuple(min(a, b) for a, b in zip(min_shape_intensities, intensity.shape))
        
        print(f"Ensuring consistent shapes: waves {min_shape_waves}, intensities {min_shape_intensities}")
        
        # Crop all arrays to minimum size
        def crop_to_shape(array, target_shape):
            """Crop array to target shape from center."""
            if array.shape == target_shape:
                return array
            
            slices = []
            for i, (current, target) in enumerate(zip(array.shape, target_shape)):
                if current > target:
                    start = (current - target) // 2
                    slices.append(slice(start, start + target))
                else:
                    slices.append(slice(None))
            
            return array[tuple(slices)]
        
        # Crop waves
        self.iam_waves = [crop_to_shape(wave, min_shape_waves) for wave in self.iam_waves]
        self.total_waves = [crop_to_shape(wave, min_shape_waves) for wave in self.total_waves]
        
        # Crop intensities
        self.iam_intensities = [crop_to_shape(intensity, min_shape_intensities) for intensity in self.iam_intensities]
        self.total_intensities = [crop_to_shape(intensity, min_shape_intensities) for intensity in self.total_intensities]
    
    def create_abtem_diffraction_objects(self):
        """
        Create abTEM diffraction objects from averaged results.
        
        Returns
        -------
        tuple
            (iam_elastic, iam_total, total_elastic, dft_total, stacked) abTEM diffraction objects
        """
        if self.iam_elastic is None or self.iam_diffraction_template is None:
            raise ValueError("No results available. Run process_all_configurations first.")
        
        # Create diffraction objects for each component
        # iam_elastic_obj = self.iam_diffraction_template.copy()
        # iam_elastic_obj.array[:] = self.iam_elastic
        
        iam_total_obj = self.iam_diffraction_template.copy()
        iam_total_obj.array[:] = self.iam_total_intensity
        
        # total_elastic_obj = self.iam_diffraction_template.copy()
        # total_elastic_obj.array[:] = self.dft_elastic
        
        dft_total_obj = self.iam_diffraction_template.copy()
        dft_total_obj.array[:] = self.dft_total_intensity
        
        # Stack for comparison
        stacked = abtem.stack([iam_total_obj,  dft_total_obj], 
                             ( "IAM Total", "DFT Total"))
        
        return iam_total_obj, dft_total_obj, stacked
    
    def visualize_diffraction_spots(self, threshold: float = 5e-5):
        """
        Visualize diffraction spots using abTEM's built-in visualization.
        
        Parameters
        ----------
        threshold : float
            Threshold for spot detection
            
        Returns
        -------
        matplotlib.figure.Figure
            The visualization figure
        """
        if self.atoms_ortho is None:
            raise ValueError("No atoms object available. Run process_all_configurations first.")
        
        # Create abTEM diffraction objects
        iam_total_obj, dft_total_obj, stacked = self.create_abtem_diffraction_objects()
        
        # Index diffraction spots
        print(f"Indexing diffraction spots (threshold = {threshold})...")
        diffraction_spots = stacked.index_diffraction_spots(self.atoms_ortho)
        
        # Use abTEM's built-in visualization
        print("Creating diffraction spot visualization...")
        visualization = diffraction_spots.show(
            explode=True, 
            common_color_scale=False, 
            cbar=True, 
            figsize=(14, 6), 
            scale=0.1
        )
        
        return visualization
    
    def show_comparison_visualization(self, threshold: float = 5e-5):
        """
        Show both abTEM diffraction spots and standard analysis side by side.
        
        Parameters
        ----------
        threshold : float
            Threshold for spot detection
        """
        print("Creating comprehensive visualization...")
        
        # Create abTEM diffraction objects
        iam_total_obj, dft_total_obj, stacked = self.create_abtem_diffraction_objects()
        
        # Show stacked diffraction patterns
        print("Showing diffraction pattern comparison...")
        stacked_viz = stacked.show(
            explode=True, 
            common_color_scale=True, 
            cbar=True, 
            figsize=(12, 6)
        )
        
        # Index and show diffraction spots
        print("Indexing and showing diffraction spots...")
        diffraction_spots = stacked.index_diffraction_spots(self.atoms_ortho)
        spots_viz = diffraction_spots.show(
            explode=True, 
            common_color_scale=False, 
            cbar=True, 
            figsize=(14, 6), 
            scale=0.1
        )
        
        # Also show standard analysis
        print("Showing phonon analysis...")
        analysis_fig = self.visualize_results(figsize=(15, 10))
        
        return stacked_viz, spots_viz, analysis_fig
    
    def visualize_elastic_inelastic_results(self, figsize: Tuple[int, int] = (18, 12)) -> plt.Figure:
        """
        Visualize elastic vs inelastic scattering separation.
        
        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        if self.iam_elastic is None:
            raise ValueError("No results to visualize. Run process_all_configurations first.")
        
        # Ensure extent is in correct format
        if self.extent is None or len(self.extent) != 4:
            ny, nx = self.iam_elastic.shape
            self.extent = [0, nx, 0, ny]
        
        fig, axes = plt.subplots(3, 3, figsize=figsize)
        
        # Row 1: IAM components
        im1 = axes[0, 0].imshow(self.iam_elastic, extent=self.extent, origin='lower', cmap='viridis')
        axes[0, 0].set_title('IAM Elastic\n|⟨ψ⟩|²')
        axes[0, 0].set_xlabel('kx (1/Å)')
        axes[0, 0].set_ylabel('ky (1/Å)')
        plt.colorbar(im1, ax=axes[0, 0])
        
        im2 = axes[0, 1].imshow(self.iam_total_intensity, extent=self.extent, origin='lower', cmap='viridis')
        axes[0, 1].set_title('IAM Total\n⟨|ψ|²⟩')
        axes[0, 1].set_xlabel('kx (1/Å)')
        axes[0, 1].set_ylabel('ky (1/Å)')
        plt.colorbar(im2, ax=axes[0, 1])
        
        im3 = axes[0, 2].imshow(self.iam_phonon, extent=self.extent, origin='lower', cmap='plasma')
        axes[0, 2].set_title('IAM Phonon\n⟨|ψ|²⟩ - |⟨ψ⟩|²')
        axes[0, 2].set_xlabel('kx (1/Å)')
        axes[0, 2].set_ylabel('ky (1/Å)')
        plt.colorbar(im3, ax=axes[0, 2])
        
        # Row 2: DFT (Total) components
        im4 = axes[1, 0].imshow(self.dft_elastic, extent=self.extent, origin='lower', cmap='viridis')
        axes[1, 0].set_title('DFT Elastic\n|⟨ψ⟩|²')
        axes[1, 0].set_xlabel('kx (1/Å)')
        axes[1, 0].set_ylabel('ky (1/Å)')
        plt.colorbar(im4, ax=axes[1, 0])
        
        im5 = axes[1, 1].imshow(self.dft_total_intensity, extent=self.extent, origin='lower', cmap='viridis')
        axes[1, 1].set_title('DFT Total\n⟨|ψ|²⟩')
        axes[1, 1].set_xlabel('kx (1/Å)')
        axes[1, 1].set_ylabel('ky (1/Å)')
        plt.colorbar(im5, ax=axes[1, 1])
        
        im6 = axes[1, 2].imshow(self.total_phonon, extent=self.extent, origin='lower', cmap='plasma')
        axes[1, 2].set_title('DFT Phonon\n⟨|ψ|²⟩ - |⟨ψ⟩|²')
        axes[1, 2].set_xlabel('kx (1/Å)')
        axes[1, 2].set_ylabel('ky (1/Å)')
        plt.colorbar(im6, ax=axes[1, 2])
        
        # Row 3: Differences (DFT - IAM)
        elastic_diff = self.dft_elastic - self.iam_elastic
        total_diff = self.dft_total_intensity - self.iam_total_intensity
        phonon_diff = self.total_phonon - self.iam_phonon
        
        diff_max1 = max(abs(elastic_diff.min()), abs(elastic_diff.max()))
        im7 = axes[2, 0].imshow(elastic_diff, extent=self.extent, origin='lower', 
                               cmap='RdBu_r', vmin=-diff_max1, vmax=diff_max1)
        axes[2, 0].set_title('Elastic Difference\n(DFT - IAM)')
        axes[2, 0].set_xlabel('kx (1/Å)')
        axes[2, 0].set_ylabel('ky (1/Å)')
        plt.colorbar(im7, ax=axes[2, 0])
        
        diff_max2 = max(abs(total_diff.min()), abs(total_diff.max()))
        im8 = axes[2, 1].imshow(total_diff, extent=self.extent, origin='lower', 
                               cmap='RdBu_r', vmin=-diff_max2, vmax=diff_max2)
        axes[2, 1].set_title('Total Difference\n(DFT - IAM)')
        axes[2, 1].set_xlabel('kx (1/Å)')
        axes[2, 1].set_ylabel('ky (1/Å)')
        plt.colorbar(im8, ax=axes[2, 1])
        
        diff_max3 = max(abs(phonon_diff.min()), abs(phonon_diff.max()))
        im9 = axes[2, 2].imshow(phonon_diff, extent=self.extent, origin='lower', 
                               cmap='RdBu_r', vmin=-diff_max3, vmax=diff_max3)
        axes[2, 2].set_title('Phonon Difference\n(DFT - IAM)')
        axes[2, 2].set_xlabel('kx (1/Å)')
        axes[2, 2].set_ylabel('ky (1/Å)')
        plt.colorbar(im9, ax=axes[2, 2])
        
        fig.suptitle(f'Elastic vs Inelastic Scattering Analysis\n'
                    f'{len(self.iam_waves)} configurations, E = {self.energy/1000:.0f} keV',
                    fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def visualize_results(self, figsize: Tuple[int, int] = (15, 10)) -> plt.Figure:
        """
        Visualize phonon-averaged diffraction results (legacy method).
        
        Returns
        -------
        matplotlib.figure.Figure
            The created figure
        """
        # Use the new elastic/inelastic visualization
        return self.visualize_elastic_inelastic_results(figsize)
    
    def get_statistics(self) -> Dict:
        """Get summary statistics for elastic/inelastic analysis."""
        if self.iam_elastic is None:
            return {}
        
        # Calculate elastic fractions
        iam_elastic_fraction = np.sum(self.iam_elastic) / np.sum(self.iam_total_intensity)
        total_elastic_fraction = np.sum(self.dft_elastic) / np.sum(self.dft_total_intensity)
        
        # Calculate phonon fractions
        iam_phonon_fraction = np.sum(self.iam_phonon) / np.sum(self.iam_total_intensity)
        total_phonon_fraction = np.sum(self.total_phonon) / np.sum(self.dft_total_intensity)
        
        # Bonding effects
        elastic_bonding_effect = np.sum(self.dft_elastic) - np.sum(self.iam_elastic)
        total_bonding_effect = np.sum(self.dft_total_intensity) - np.sum(self.iam_total_intensity)
        phonon_bonding_effect = np.sum(self.total_phonon) - np.sum(self.iam_phonon)
        
        return {
            'n_configs': len(self.iam_waves),
            
            # Total intensities
            'iam_total_intensity': np.sum(self.iam_total_intensity),
            'dft_total_intensity': np.sum(self.dft_total_intensity),
            
            # Elastic components
            'iam_elastic_intensity': np.sum(self.iam_elastic),
            'dft_elastic_intensity': np.sum(self.dft_elastic),
            'iam_elastic_fraction': iam_elastic_fraction,
            'dft_elastic_fraction': total_elastic_fraction,
            
            # Phonon components
            'iam_phonon_intensity': np.sum(self.iam_phonon),
            'dft_phonon_intensity': np.sum(self.total_phonon),
            'iam_phonon_fraction': iam_phonon_fraction,
            'dft_phonon_fraction': total_phonon_fraction,
            
            # Bonding effects
            'elastic_bonding_effect': elastic_bonding_effect,
            'total_bonding_effect': total_bonding_effect,
            'phonon_bonding_effect': phonon_bonding_effect,
            
            # RMS measures
            'elastic_diff_rms': np.sqrt(np.mean((self.dft_elastic - self.iam_elastic)**2)),
            'phonon_diff_rms': np.sqrt(np.mean((self.total_phonon - self.iam_phonon)**2)),
            'total_diff_rms': np.sqrt(np.mean((self.dft_total_intensity - self.iam_total_intensity)**2))
        }
    
    def save_results(self, output_dir: str = "phonon_results"):
        """Save results to files."""
        import os
        import json
        
        os.makedirs(output_dir, exist_ok=True)
        
        if self.iam_elastic is not None:
            # Save elastic components
            np.save(f"{output_dir}/iam_elastic.npy", self.iam_elastic)
            np.save(f"{output_dir}/dft_elastic.npy", self.dft_elastic)
            
            # Save total intensities
            np.save(f"{output_dir}/iam_total.npy", self.iam_total_intensity)
            np.save(f"{output_dir}/dft_total.npy", self.dft_total_intensity)
            
            # Save phonon components
            np.save(f"{output_dir}/iam_phonon.npy", self.iam_phonon)
            np.save(f"{output_dir}/dft_phonon.npy", self.total_phonon)
            
            # Save metadata
            metadata = {
                'n_configs': len(self.iam_waves),
                'sampling': self.sampling.tolist() if self.sampling is not None else None,
                'extent': self.extent,
                'energy': self.energy,
                'max_angle': self.max_angle,
                'statistics': self.get_statistics()
            }
            
            with open(f"{output_dir}/metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2, default=str)
        
        print(f"Results saved to {output_dir}")