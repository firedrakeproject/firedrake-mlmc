all:

lint:
	@echo "    Linting firedrake-mlmc codebase"
	@flake8 firedrake_mlmc
	@echo "    Linting firedrake-mlmc test suite"
	@flake8 tests
	@echo "    Linting firedrake-mlmc demo suite"
	@flake8 examples

test:
	@echo "    Running all tests"
	@py.test tests $(PYTEST_ARGS)
