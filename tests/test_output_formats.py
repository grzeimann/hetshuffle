from __future__ import annotations

import io

from hetdex_shuffle.obj_dist_ifu_vis import writeFailed, writeHeader, writeShuffle


def test_distance_output_header_and_success_row():
    output = io.StringIO()
    writeHeader(output)
    writeShuffle(output, 180.0, 30.0, 1, [58], 0.25, -0.5,
                 2, [66], 0.75, -1.0, 0)

    text = output.getvalue()
    assert "#RA" in text
    assert "180.000000" in text
    assert "0.25" in text


def test_distance_output_failure_row_is_descriptive():
    output = io.StringIO()
    writeFailed(output, 180.0, 30.0, -1)

    text = output.getvalue()
    assert "180.000000" in text
    assert "No IFUs within close proximity." in text
