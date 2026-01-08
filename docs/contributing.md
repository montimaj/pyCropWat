# Contributing

Thank you for your interest in contributing to pyCropWat!

## Development Setup

1. Fork and clone the repository:
   ```bash
   git clone https://github.com/montimaj/pyCropWat.git
   cd pyCropWat
   ```

2. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

4. Authenticate with Google Earth Engine:
   ```bash
   earthengine authenticate
   ```

## Code Style

We use the following tools to maintain code quality:

- **Black** for code formatting
- **Ruff** for linting
- **NumPy-style** docstrings

Run formatting and linting:

```bash
# Format code
black pycropwat tests

# Lint code
ruff check pycropwat tests
```

## Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pycropwat --cov-report=html

# Run specific test file
pytest tests/test_core.py
```

## Documentation

Documentation is built with MkDocs. To preview locally:

```bash
pip install mkdocs mkdocs-material mkdocstrings[python]
mkdocs serve
```

Then open http://localhost:8000 in your browser.

## Pull Request Process

1. Create a new branch for your feature:
   ```bash
   git checkout -b feature/my-new-feature
   ```

2. Make your changes and commit:
   ```bash
   git add .
   git commit -m "Add my new feature"
   ```

3. Push to your fork:
   ```bash
   git push origin feature/my-new-feature
   ```

4. Open a Pull Request on GitHub

### PR Checklist

- [ ] Code follows the project style guidelines
- [ ] Tests pass locally
- [ ] New functionality includes tests
- [ ] Documentation updated if needed
- [ ] Commit messages are clear and descriptive

## Reporting Issues

When reporting issues, please include:

- Python version
- pyCropWat version
- Operating system
- Full error traceback
- Minimal code to reproduce the issue

## Feature Requests

Feature requests are welcome! Please:

- Check existing issues first
- Describe the use case clearly
- Explain why this would benefit other users

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
