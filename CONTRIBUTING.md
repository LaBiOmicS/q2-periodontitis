# Contributing to `q2-periodontitis`

Thank you for your interest in contributing to `q2-periodontitis`! We welcome contributions from bioinformaticians, oral microbiologists, developers, and researchers.

---

## 🚀 How to Contribute

### 1. Reporting Bugs
- Search existing [GitHub Issues](https://github.com/LaBiOmicS/q2-periodontitis/issues) to verify the bug hasn't been reported.
- Open a new issue using the **Bug Report** template.
- Include full CLI command invocations, QIIME 2 version, and complete error tracebacks.

### 2. Suggesting Features
- Open an issue using the **Feature Request** template.
- Provide biological context, mathematical rationale, or clinical use cases for the requested feature.

### 3. Submitting Pull Requests (PRs)
1. Fork the repository and create a feature branch:
   ```bash
   git checkout -b feature/my-new-feature
   ```
2. Install the package in editable mode with development dependencies:
   ```bash
   pip install -e .
   ```
3. Write clean, PEP 8 compliant code and add unit tests under `q2_periodontitis/tests/`.
4. Run unit tests locally using `pytest`:
   ```bash
   pytest
   ```
5. Commit your changes with informative commit messages and push to your fork.
6. Submit a Pull Request targeting the `main` branch.

---

## 📜 Code Style & Testing
- Code formatting follows PEP 8 standards.
- Ensure all public functions have descriptive docstrings.
- PRs will only be merged if all GitHub Actions CI checks pass.
