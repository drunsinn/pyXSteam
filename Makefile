.PHONY: clean

clean:
	rm -rf build dist *.egg-info

test:
	uv run pytest

spell:
	uv run codespell .

lint:
	uv run ruff check

format:
#	uv run ruff format
	uv run black . --config ./pyproject.toml

build:
	uv build

doc:
	cd docs && $(MAKE) html

all: clean test spell lint format build doc



