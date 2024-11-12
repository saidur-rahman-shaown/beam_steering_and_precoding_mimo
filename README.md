
# Precoding and Beam Steering in a MIMO-OFDM System Using Ray Tracing in Sub-6GHz

This project explores beam steering techniques in a MIMO-OFDM system within the sub-6GHz frequency band using ray tracing simulations. The setup features a transmitter placed Electrical and Computer Engineering building and a receiver located several meters away. The simulation focuses on evaluating bit error rates (BER) and optimizing beam alignment for effective communication.

## Features
- **MIMO-OFDM Precoding**: Models bit error rates for a multiple-input multiple-output orthogonal frequency division multiplexing system in the sub-6GHz spectrum.
- **Ray Tracing Integration**: Implements ray tracing techniques to achieve precise beam steering and enhance signal alignment.
- **Realistic Setup**: Simulates practical deployment scenarios for wireless communication.
- **Customizable Parameters**: Easily adapt system configurations for different environments and use cases.

## Prerequisites
- MATLAB (R2020a or later is recommended)
- Communications Toolbox
- Phased Array System Toolbox
- Anternna toolbox

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/Beam-Steering-MIMO-OFDM.git
   ```
2. Open MATLAB and navigate to the project folder.
3. Ensure the required toolboxes are installed.

## Usage
1. Configure array parameters in the `array_setup.m` file:
   - Transmitter and receiver positions

2. Configure system parameters in `MIMOOfdmBeamSteering.m` file:
   - Number of antennas and subcarriers
   - Frequency and bandwidth settings
2. Run the main simulation script:
   ```matlab
   run('main_simulation.m')
   ```
3. Analyze the results:
   - BER vs. SNR graphs
   - Beamforming patterns
   - Log files detailing system performance

## Project Structure
- **main_simulation.m**: Main script for running the simulation.
- **config.m**: Contains configurable parameters for the simulation.
- **raytracing/**: Ray tracing utilities for path calculations.
- **mimo_ofdm/**: Core MIMO-OFDM algorithms and helper functions.
- **results/**: Stores simulation outputs such as plots and logs.

## Results
This simulation provides:
- BER performance curves for varying SNR levels.
- Visualizations of beamforming patterns.
- Logs detailing system configurations and outcomes.

## Future Work
- Extension to mmWave frequencies.
- Inclusion of NLOS (Non-Line-of-Sight) scenarios.
- Implementation of dynamic beam tracking for moving receivers.

## Contributing
Contributions are welcome! Here's how you can contribute:
1. Fork the repository.
2. Create a new branch for your feature (`feature/your-feature-name`).
3. Push your changes and open a pull request.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact
For any queries or feedback, please contact:
- **Author**: [Your Name]
- **Email**: your.email@example.com
