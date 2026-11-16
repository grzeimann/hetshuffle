.PHONY: help lint test coverage ci clean

help:
	@echo "Common targets:"
	@echo "  make lint      - Run Ruff lint checks"
	@echo "  make test      - Run the pytest suite"
	@echo "  make coverage  - Run pytest with coverage reports"
	@echo "  make ci        - Run the checks used by CI"
	@echo "  make clean     - Remove test and coverage artifacts"

lint:
	python -m ruff check tests

test:
	python -m pytest -q

coverage:
	python -m coverage run -m pytest -q
	python -m coverage report -m
	python -m coverage xml -o coverage.xml

ci: lint coverage

clean:
	rm -rf .pytest_cache .coverage coverage.xml htmlcov
