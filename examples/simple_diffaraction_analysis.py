#!/usr/bin/env python3
"""
Simple Phonon Analysis
======================
script for phonon-averaged diffraction analysis for IAM and IAM+bonding.
"""

import matplotlib.pyplot as plt
from src.phonon_diffraction import PhononDiffraction


def main():
    """Main analysis function."""
    
    print("="*50)
    print("ELASTIC vs INELASTIC SCATTERING ANALYSIS")
    print("="*50)
    
    # Initialize analyzer
    analyzer = PhononDiffraction(
        data_dir="pot-diff/pot_diff/data",
        energy=100e3,  # 100 keV
        max_angle=50   # 50 mrad
    )
    
    # Process all configurations
    summary = analyzer.process_all_configurations()
    
    if summary['successful'] == 0:
        print("No configurations processed successfully!")
        return
    
    # Get statistics
    stats = analyzer.get_statistics()
    
    # Print results (diffraction patterns only)
    print(f"\nResults:")
    print(f"  Configurations: {stats['n_configs']}")
    print(f"  IAM total intensity: {stats['iam_total_intensity']:.2e}")
    print(f"  DFT total intensity: {stats['dft_total_intensity']:.2e}")
    print(f"  Total bonding effect: {stats['total_bonding_effect']:.2e}")
    print(f"  Total RMS difference: {stats['total_diff_rms']:.2e}")
    print(f"  Elastic RMS difference: {stats['elastic_diff_rms']:.2e}")
    print(f"  Phonon RMS difference: {stats['phonon_diff_rms']:.2e}")
    
    # Visualize diffraction spots only
    print("\nCreating diffraction spot visualization...")
    try:
        visualization = analyzer.visualize_diffraction_spots(threshold=5e-5)
        print("✅ Diffraction spots visualization created")
        
        # Save the diffraction spots visualization
        analyzer.save_results("simple_phonon_results")
        if hasattr(visualization, 'savefig'):
            visualization.savefig("simple_phonon_results/diffraction_spots.png", dpi=150, bbox_inches='tight')
        
    except Exception as e:
        print(f"Error: Could not create diffraction spots visualization: {e}")
        print("This might be due to shape mismatch or missing abTEM objects.")
        return
    
    print(f"\nAnalysis complete!")
    print(f"Results saved to: simple_phonon_results/")
    print(f"Diffraction spots visualization: simple_phonon_results/diffraction_spots.png")
    print(f"Diffraction spots visualization displayed above")
    
    plt.show()


def test_single_config():
    """Test with single configuration."""
    print("Testing with single configuration...")
    
    analyzer = PhononDiffraction()
    summary = analyzer.process_all_configurations(max_configs=1)
    
    if summary['successful'] > 0:
        fig = analyzer.visualize_results()
        fig.savefig("test_single_config.png", dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("✓ Single config test successful!")
        return True
    else:
        print("✗ Single config test failed!")
        return False


if __name__ == "__main__":
    # Uncomment to test with single config first
    # test_single_config()
    
    # Run full analysis
    main()