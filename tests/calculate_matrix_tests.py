import os
import subprocess
import tempfile
from pathlib import Path

import pytest
from hmmcdr.cli import \
    calculate_matrices  # Adjust import based on your package structure


# Create test data fixtures
@pytest.fixture
def sample_input_file():
    content = """chr17_PAT\t1000\t2000\t0.8
chr17_PAT\t2000\t3000\t0.9
chr17_PAT\t3000\t4000\t0.7"""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
        f.write(content)
    yield f.name
    os.unlink(f.name)


@pytest.fixture
def reference_file():
    content = """chr17_PAT\t1000\t2000\tregion1
chr17_PAT\t2000\t3000\tregion2
chr17_PAT\t3000\t4000\tregion3"""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as f:
        f.write(content)
    yield f.name
    os.unlink(f.name)


@pytest.fixture
def output_file():
    with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as f:
        yield f.name
        os.unlink(f.name)


def test_calculate_matrices_cli(sample_input_file, reference_file, output_file):
    """Test the CLI interface of calculate_matrices"""
    cmd = ["calculate_matrices", sample_input_file, reference_file, output_file]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # Verify output file exists
    assert os.path.exists(output_file)

    # Verify output file content
    with open(output_file) as f:
        content = f.read()
        # Add assertions based on expected output format
        assert len(content.strip().split("\n")) > 0


def test_calculate_matrices_function(sample_input_file, reference_file, output_file):
    """Test the Python function directly"""
    result = calculate_matrices(
        input_file=sample_input_file,
        reference_file=reference_file,
        output_file=output_file,
    )

    assert os.path.exists(output_file)
    # Add more specific assertions based on expected output


def test_invalid_input_file():
    """Test behavior with non-existent input file"""
    with pytest.raises(FileNotFoundError):
        calculate_matrices(
            input_file="nonexistent.bed",
            reference_file="reference.bed",
            output_file="output.bed",
        )


def test_invalid_file_format(tmp_path):
    """Test behavior with malformed input file"""
    invalid_file = tmp_path / "invalid.bed"
    invalid_file.write_text("invalid content\nmore invalid content")

    with pytest.raises(ValueError):
        calculate_matrices(
            input_file=str(invalid_file),
            reference_file="reference.bed",
            output_file="output.bed",
        )


def test_output_content(sample_input_file, reference_file, output_file):
    """Test the actual content of the output file"""
    calculate_matrices(
        input_file=sample_input_file,
        reference_file=reference_file,
        output_file=output_file,
    )

    with open(output_file) as f:
        content = f.readlines()

    # Add assertions specific to your output format
    # For example:
    # assert content[0].startswith("chr17_PAT")
    # assert len(content) == expected_number_of_lines
    # assert all required columns are present
    # assert values are within expected ranges
