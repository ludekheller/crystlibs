# Contributing to Crystallography Python Libraries

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## 📋 Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)

## 🤝 Code of Conduct

This project adheres to a code of conduct that all contributors are expected to follow:

- Be respectful and inclusive
- Accept constructive criticism gracefully
- Focus on what is best for the community
- Show empathy towards other community members

## 🚀 Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/crystallography-python.git
   cd crystallography-python
   ```
3. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## 💡 How to Contribute

### Reporting Bugs

If you find a bug, please create an issue with:
- Clear, descriptive title
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Your environment (Python version, OS, library versions)
- Minimal code example that demonstrates the issue

### Suggesting Enhancements

For feature requests:
- Use a clear, descriptive title
- Provide detailed description of the proposed feature
- Explain why this enhancement would be useful
- Include code examples if applicable

### Pull Requests

1. Ensure your code follows the project's coding standards
2. Add tests for new functionality
3. Update documentation as needed
4. Ensure all tests pass
5. Submit your pull request

## 🔧 Development Setup

### Install Development Dependencies

```bash
pip install -r requirements.txt
pip install pytest pytest-cov flake8 black
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=. --cov-report=html

# Run specific test file
pytest tests/test_crystlib.py
```

### Code Formatting

We use `black` for code formatting:

```bash
# Format all Python files
black *.py

# Check formatting without making changes
black --check *.py
```

### Linting

We use `flake8` for code linting:

```bash
# Lint all files
flake8 *.py

# Lint specific file
flake8 crystlib.py
```

## 📝 Coding Standards

### General Guidelines

- **Python Version**: Code should be compatible with Python 2.7 and 3.6+
- **PEP 8**: Follow PEP 8 style guide (enforced by flake8)
- **Line Length**: Maximum 88 characters (black default)
- **Docstrings**: All functions must have comprehensive docstrings

### Docstring Format

```python
def function_name(param1, param2):
    """
    One-line description of the function.
    
    Longer description explaining the function's purpose,
    algorithm, and any important notes.
    
    Input:
        param1: type - Description of param1
        param2: type - Description of param2
    
    Output:
        return_type - Description of return value
    
    Usage Example:
        >>> import crystlib
        >>> result = crystlib.function_name(arg1, arg2)
        >>> print(result)
        Expected output
    
    Notes:
        - Important note 1
        - Important note 2
    
    Formula:
        Mathematical formula if applicable
    """
    # Function implementation
    pass
```

### Naming Conventions

- **Functions**: Use lowercase with underscores: `calculate_distance()`
- **Variables**: Use lowercase with underscores: `lattice_vector`
- **Constants**: Use uppercase with underscores: `MAX_ITERATIONS`
- **Classes**: Use CamelCase: `LatticePlotter`

### NumPy Functions

- Prefix NumPy-optimized versions with `np_`: `np_euler_matrix()`
- Keep list-based versions for compatibility: `euler_matrix()`

## 🧪 Testing

### Test Organization

Tests are organized in the `tests/` directory:
- `test_crystlib.py` - Tests for crystlib functions
- `test_orilib.py` - Tests for orilib functions
- `test_plotlib.py` - Tests for plotlib functions
- `test_projlib.py` - Tests for projlib functions

### Writing Tests

```python
import pytest
import numpy as np
from crystlib import cubic_lattice_vec

def test_cubic_lattice_vec():
    """Test cubic lattice vector generation."""
    a = 3.6
    lattice = cubic_lattice_vec(a)
    
    # Check dimensions
    assert lattice.shape == (3, 3)
    
    # Check values
    expected = np.diag([a, a, a])
    np.testing.assert_array_almost_equal(lattice, expected)
    
def test_cubic_lattice_vec_negative():
    """Test that negative lattice parameter raises error."""
    with pytest.raises(ValueError):
        cubic_lattice_vec(-3.6)
```

### Test Coverage

- Aim for >80% code coverage
- Test both normal and edge cases
- Test error handling
- Include numerical precision tests

## 📖 Documentation

### Updating Documentation

When adding new functions:
1. Add comprehensive docstring to the function
2. Update the appropriate `*_DOCUMENTATION_COMPLETE.md` file
3. Update the `*_QUICK_REFERENCE.md` file
4. Add usage examples

### Documentation Structure

Each library has 4 documentation levels:
- **COMPLETE**: Full reference with examples
- **SUMMARY**: Organized overview
- **QUICK_REFERENCE_COMPREHENSIVE**: Detailed quick reference
- **QUICK_REFERENCE**: Ultra-fast reference

### Adding Examples

Place example scripts in the `examples/` directory:
```python
# examples/my_example.py
"""
Example: Demonstrating new feature

This example shows how to use the new functionality.
"""
import crystlib

# Your example code here
```

## 📤 Submitting Changes

### Before Submitting

- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] New tests added for new functionality
- [ ] Documentation updated
- [ ] Code formatted with black
- [ ] No linting errors
- [ ] Commit messages are clear and descriptive

### Commit Messages

Use clear, descriptive commit messages:

```
Add function for hexagonal slip system generation

- Implement gensystemsHex() function
- Add comprehensive docstring
- Include usage examples
- Add tests for edge cases
```

### Pull Request Process

1. **Update your branch** with the latest main:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Push your changes**:
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Create Pull Request** on GitHub with:
   - Clear title describing the change
   - Detailed description of what and why
   - Reference any related issues
   - List of changes made
   - Screenshots if applicable

4. **Address review comments**:
   - Respond to all comments
   - Make requested changes
   - Push updates to same branch

5. **Merge**: Once approved, maintainers will merge your PR

## 🎯 Areas for Contribution

### High Priority

- [ ] Additional unit tests for edge cases
- [ ] Performance optimization for large datasets
- [ ] Additional crystal system support
- [ ] Integration with common data formats (HDF5, DREAM.3D)

### Medium Priority

- [ ] Interactive Jupyter notebook tutorials
- [ ] Additional visualization options
- [ ] More comprehensive examples
- [ ] Performance benchmarking

### Good First Issues

- [ ] Documentation improvements
- [ ] Adding type hints
- [ ] Improving error messages
- [ ] Adding docstring examples

## 📚 Resources

- [NumPy Documentation](https://numpy.org/doc/)
- [SciPy Documentation](https://docs.scipy.org/)
- [Matplotlib Documentation](https://matplotlib.org/stable/contents.html)
- [Python Testing with pytest](https://docs.pytest.org/)
- [PEP 8 Style Guide](https://www.python.org/dev/peps/pep-0008/)

## ❓ Questions?

If you have questions about contributing:
- Open a [GitHub Discussion](https://github.com/yourusername/crystallography-python/discussions)
- Create an issue with the "question" label
- Contact the maintainers

## 🙏 Thank You!

Your contributions help make this project better for the entire materials science community!

---

*Last Updated: December 2024*
