
# Flame Components

This repository provides several functions for calculating various fire behavior metrics such as mid-flame wind speed, flame length, flame height, flame tilt, flame residence time, and flame depth. It also includes utility functions for multiprocessing these calculations across blocks of data, which can be useful for high-performance fire modeling.

## Features

- **Mid-Flame Wind Speed Calculation** (`getMidFlameWS`): Calculates mid-flame wind speed based on parameters such as wind speed, canopy cover, and canopy height.
- **Flame Length Estimation** (`getFlameLength`): Estimates flame length using different published models.
- **Flame Height Calculation** (`getFlameHeight`): Calculates flame height based on flame length and model-specific parameters.
- **Flame Tilt Angle Calculation** (`getFlameTilt`): Calculates the angle of flame tilt relative to vertical.
- **Flame Residence Time** (`getFlameResidenceTime`): Computes flame residence time based on rate of spread, fuel consumption, and wind speed.
- **Flame Depth Calculation** (`getFlameDepth`): Calculates flame depth using flame residence time and rate of spread.
- **Array Multiprocessing** (`flameComponent_ArrayMultiprocessing`): Enables multiprocessing of flame component calculations across blocks of data, making large-scale processing efficient.

## Requirements

- Python 3.8+
- **Libraries**: `numpy`, `multiprocessing`

## Usage

### Key Functions

- **`getMidFlameWS`**: Calculates mid-flame wind speed.
- **`getFlameLength`**: Estimates flame length based on specified models.
- **`getFlameHeight`**: Calculates flame height for a given flame length.
- **`getFlameTilt`**: Computes flame tilt angle using various models.
- **`getFlameResidenceTime`**: Estimates flame residence time for a given rate of spread.
- **`getFlameDepth`**: Computes flame depth from flame residence time and rate of spread.

### Example

```python
from flame_components import getMidFlameWS, getFlameLength

# Example calculation for mid-flame wind speed
mid_flame_ws = getMidFlameWS(
    wind_speed=15,
    canopy_cover=50,
    canopy_ht=10,
    canopy_baseht=2,
    units='SI'
)

# Example calculation for flame length
flame_length = getFlameLength(
    model='Byram_HEAD',
    fire_intensity=500
)
```

### Array Multiprocessing Example

To perform calculations across blocks of data using multiple processors:

```python
from flame_components import flameComponent_ArrayMultiprocessing

# Example multiprocessing calculation
results = flameComponent_ArrayMultiprocessing(
    flame_function='midflame_ws',
    num_processors=4,
    wind_speed=array_of_wind_speed,
    canopy_cover=array_of_canopy_cover,
    canopy_ht=array_of_canopy_ht,
    canopy_baseht=array_of_canopy_baseht,
    units='SI'
)
```

## Contributing

Contributions are welcome. If you have ideas for new features or improvements, please submit a pull request or open an issue.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
