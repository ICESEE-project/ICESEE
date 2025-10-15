# Zarr Setup Tests

Test suite for Zarr storage configuration and functionality in ICESEE.

## Purpose

Validate Zarr-based storage for efficient, cloud-optimized ensemble data management.

## Overview

Zarr provides:
- Chunked array storage
- Compression options
- Cloud storage compatibility
- Parallel I/O support
- Efficient partial reads/writes

These tests ensure ICESEE's Zarr integration works correctly.

## Test Coverage

Tests include:
- Zarr array creation
- Data writing and reading
- Chunk size optimization
- Compression settings
- Metadata handling
- Parallel I/O with Zarr
- Cloud storage integration (if configured)

## Running Tests

### Local Tests
```bash
cd src/tests/zarr_setup
python -m pytest
```

### With Specific Storage Backend
```bash
# Test with local filesystem
python -m pytest --storage local

# Test with S3 (if configured)
python -m pytest --storage s3
```

## Configuration

Zarr settings can be configured in test fixtures:
```python
zarr_config = {
    'chunks': (10, 100, 100),
    'compressor': 'blosc',
    'compression_level': 5,
}
```

## Performance Tests

### Chunk Size Optimization
```bash
python benchmark_chunk_sizes.py
```

### Compression Comparison
```bash
python test_compression_methods.py --output results.csv
```

## Storage Backends

Tests support multiple backends:
- **Local filesystem**: Default for development
- **S3**: AWS S3 or compatible (MinIO, etc.)
- **Network filesystems**: Lustre, GPFS

## Common Operations

### Creating Zarr Store
```python
import zarr
store = zarr.DirectoryStore('data.zarr')
root = zarr.group(store=store)
ensemble = root.create_dataset('ensemble', shape=(50, 1000), chunks=(5, 100))
```

### Reading Data
```python
ensemble = zarr.open('data.zarr/ensemble', mode='r')
member_0 = ensemble[0, :]
```

## Best Practices

- Choose chunk sizes based on access patterns
- Use compression for I/O-bound workloads
- Test with realistic data sizes
- Profile I/O performance
- Consider cloud storage for large datasets

## Troubleshooting

### Slow I/O
- Check chunk sizes
- Verify compression settings
- Monitor network/disk bandwidth

### Compatibility Issues
- Ensure Zarr version compatibility
- Check storage backend support

## Requirements

- zarr (Python package)
- numcodecs (compression)
- fsspec (filesystem interfaces)
- Optional: s3fs (for S3 support)

## Expected Results

Tests validate:
- Data integrity after write/read cycles
- Correct handling of array shapes
- Proper compression/decompression
- Efficient chunked access
- Metadata preservation

## Adding Tests

When adding Zarr tests:
1. Test both read and write operations
2. Verify data integrity
3. Include performance benchmarks
4. Test edge cases (empty arrays, large arrays)
5. Clean up test data files

## Future Enhancements

- Distributed Zarr with Dask
- Cloud storage optimization
- Advanced compression codecs
- Automatic chunk size selection
