# Run Model Data Assimilation

This directory contains the main entry points and orchestration logic for running data assimilation workflows with ICESEE.

## Main Scripts

### run_models_da.py
Primary entry point for running data assimilation experiments.

**Features**:
- Configuration loading
- Workflow orchestration
- Model-EnKF integration
- Results management

**Usage**:
```bash
python run_models_da.py -F params.yaml --Nens 50
```

### icesee_da_serial.py
Serial (single-processor) data assimilation implementation.

**When to use**:
- Small problems
- Debugging
- Development
- Testing

### icesee_da_partial_parallel.py
Partially parallelized implementation.

**Features**:
- Parallel forecast step
- Serial analysis step
- Suitable for medium-scale problems

### icesee_da_full_parallel.py
Fully parallelized data assimilation implementation.

**Features**:
- Parallel forecast and analysis
- Distributed I/O
- Maximum scalability
- For production HPC runs

## Supporting Modules

### _error_generation.py
Functions for generating model and observation errors.

**Capabilities**:
- Gaussian random fields
- Spatially correlated errors
- Temporal correlation
- Observation error simulation

### _localization_functions.py
Covariance localization implementations.

**Methods**:
- Gaspari-Cohn correlation function
- Distance-based localization
- Custom localization patterns
- Adaptive localization

## Workflow Structure

A typical ICESEE workflow:

1. **Initialization**
   - Load configuration
   - Initialize model
   - Generate initial ensemble

2. **Forecast Step**
   - Run ensemble members forward in time
   - (Parallel execution)

3. **Analysis Step**
   - Compute ensemble statistics
   - Apply EnKF update
   - (Parallel or serial)

4. **Output**
   - Save ensemble state
   - Checkpoint data
   - Diagnostic output

5. **Repeat** for next observation time

## Execution Modes

### Serial
```bash
python icesee_da_serial.py --config params.yaml
```

### Partial Parallel
```bash
mpirun -np 24 python icesee_da_partial_parallel.py --config params.yaml
```

### Full Parallel
```bash
mpirun -np 48 python icesee_da_full_parallel.py --config params.yaml
```

## Configuration

Key configuration parameters (see `config/README.md` for full details):
- `Nens`: Ensemble size
- `execution_mode`: 0=serial, 1=partial parallel, 2=full parallel
- `model_name`: Which model to use
- `parallel_flag`: Enable/disable parallelization

## Error Handling

The framework includes robust error handling:
- Configuration validation
- Model execution errors
- I/O failures
- MPI communication errors

## Restart Capability

Workflows can be restarted from checkpoints:
- Automatic checkpoint creation
- Configurable checkpoint frequency
- Resume from any saved state

## Logging

Comprehensive logging for debugging and monitoring:
- Configuration summary
- Timing information
- Performance metrics
- Error messages

## Testing

Test these modules with the examples in `applications/` directories.

## Performance Tips

- Use full parallel mode for production runs
- Profile to identify bottlenecks
- Balance ensemble size with available processors
- Use localization for large spatial domains
- Enable checkpointing for long runs

## References

See the [ICESEE Wiki](https://github.com/ICESEE-project/ICESEE/wiki) for detailed workflow documentation.
