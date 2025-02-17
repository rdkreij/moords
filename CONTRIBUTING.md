# Contributing to `moords`

Thank you for considering contributing to the `moords` project! Contributions are welcome and appreciated. 

## How to Contribute

### Setting Up the Development Environment

Before you start contributing, you need to clone the repository and set up a virtual environment.

1. Clone the repository:
   ```bash
   git clone https://github.com/rdkreij/moords.git
   cd moords
   ```

2. Install [Poetry](https://python-poetry.org/) if you haven't already.
    ```bash
    pip install poetry
    ```

3. Install the dependencies using Poetry:
   ```bash
   poetry install
   ```

4. To activate the virtual environment, run:
   ```bash
   poetry shell
   ```

### Reporting Bugs

If you encounter a bug, please check if the issue has already been reported by searching the [issues page](https://github.com/rdkreij/moords/issues). If not, create a new issue with the following information:

- A clear description of the bug
- Steps to reproduce the bug
- Expected behavior
- Any relevant logs or error messages

### Suggesting Features

If you have an idea for a new feature or enhancement, feel free to create a new issue with a description of the idea. Be sure to explain:

- What the feature does
- Why it would be useful
- Any possible implementation suggestions

### Pull Requests

To submit a pull request (PR), follow these steps:

1. Fork the repository and create your branch (`git checkout -b feature-branch`).
2. Make your changes, ensuring they are well-documented and tested.
3. If possible, update relevant documentation and/or examples.
4. Push your changes to your fork (`git push origin feature-branch`).
5. Open a pull request to the main repository.

Please make sure your code follows the existing style and includes relevant tests.

### Testing Changes

We use `unittest` for testing. Please ensure that any new code or features are accompanied by appropriate tests. You can run the tests locally using:

```bash
python -m unittest discover
```

### Code Style

Please adhere to the following coding style conventions:

- **PEP 8**: Follow the Python style guide for code formatting.
- **Docstrings**: Use docstrings for all public modules, functions, and classes.
- **Type hints**: Use type hints for function signatures where applicable.

The project uses `pre-commit` hooks to ensure consistent formatting and code quality. To set this up, run the following commands to install the hooks:

```bash
pre-commit install
```

This will automatically check your code for style issues (including PEP 8) and other common problems when committing changes.

## License

By contributing to `moords`, you agree that your contributions will be licensed under the [MIT License](LICENSE).
