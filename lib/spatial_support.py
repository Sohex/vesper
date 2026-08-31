"""Versioned spatial-support contracts and their fail-closed validator.

The contract is metadata, not a remapper. ``lib.gridding`` and ``lib.remap``
own operators; this module says what support and semantics an artifact claims
those operators acted on. SPAT-1.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
DECLARATION = ROOT / "config" / "spatial_support.yaml"
CONVERSION_DECLARATION = ROOT / "config" / "spatial_conversion.yaml"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class SpatialContractError(ValueError):
    """A spatial contract is incomplete or internally contradictory."""


def declaration(path: Path = DECLARATION) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def conversion_declaration(path: Path = CONVERSION_DECLARATION) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def canonical_digest(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def grid_support_contract(spec: Any, artifact_id: str, support_id: str,
                          kind: str, source_artifacts: list[str],
                          native_measures: list[dict[str, str]] | None = None) -> dict:
    """Construct the minimal static contract for a bounded structured grid.

    ``spec`` is deliberately duck-typed to ``gridding.GridSpec`` so this
    metadata module does not import the operator module it describes.
    """
    edges = {"longitude_edges": list(map(float, spec.lon_edges)),
             "sine_latitude_edges": list(map(float, spec.sin_edges))}
    measures = native_measures or [
        {"kind": "area", "variable": "cell_area", "units": "m2"},
    ]
    contract = {
        "contract_version": declaration()["contract_version"],
        "artifact_id": artifact_id,
        "supports": [{
            "id": support_id,
            "kind": kind,
            "dimensions": {"latitude": int(spec.nlat), "longitude": int(spec.nlon)},
            "geometry_sha256": canonical_digest({
                "kind": kind, "shape": [int(spec.nlat), int(spec.nlon)],
                "constructor": str(spec.source),
            }),
            "coordinates_sha256": canonical_digest(edges),
            "native_measures": measures,
        }],
        "fields": [{
            "name": "cell_coverage",
            "semantics": "dimensionless_fraction",
            "units": "1",
            "support_id": support_id,
            "spatial_dimensions": ["latitude", "longitude"],
            "effective_measure": {"mode": "native", "measure_kind": "area"},
            "aggregation": {
                "operator": "identity",
                "operator_version": "lib.spatial_support.grid_support_contract/1",
                "source_support": support_id,
                "destination_support": support_id,
                "nonlinear_order": "not_applicable",
            },
        }],
        "time_support": {"semantics": "static"},
        "provenance": {
            "source_artifacts": source_artifacts,
            "operator_versions": {
                "identity": "lib.spatial_support.grid_support_contract/1",
            },
        },
    }
    validate_contract(contract)
    return contract


def _quantile_probabilities(field: dict, where: str) -> list[float]:
    """The probability vector a quantile table is unreadable without.

    `lib/gridding.py:cell_quantiles` requires `probs` from its caller and
    declines to default it, because a quantile vector states which part of the
    distribution a consumer needs resolved and on this world the tails are where
    the ice and the abyssal floor are. The same argument applies one layer out:
    a table whose probabilities are not written down beside it is a block of
    numbers in the field's units and nothing more, so the contract carries them
    under the same precondition the operator enforces.
    """
    probs = field.get("quantile_probabilities")
    if not isinstance(probs, list) or not probs:
        raise SpatialContractError(
            f"{where} needs quantile_probabilities; a quantile table without the "
            "probabilities it is taken at cannot be read")
    if any(not isinstance(p, (int, float)) or isinstance(p, bool) for p in probs):
        raise SpatialContractError(f"{where}.quantile_probabilities must be numbers")
    values = [float(p) for p in probs]
    if values[0] < 0.0 or values[-1] > 1.0 or \
            any(b <= a for a, b in zip(values, values[1:])):
        raise SpatialContractError(
            f"{where}.quantile_probabilities must strictly increase inside [0, 1]; "
            "an unordered vector describes an unordered table and nothing says so")
    return values


def _need(mapping: dict, key: str, where: str) -> Any:
    if key not in mapping:
        raise SpatialContractError(f"{where} is missing {key!r}")
    return mapping[key]


def _positive_dimensions(value: Any, where: str) -> dict[str, int]:
    if not isinstance(value, dict) or not value:
        raise SpatialContractError(f"{where} dimensions must be a non-empty mapping")
    if any(not isinstance(v, int) or isinstance(v, bool) or v <= 0
           for v in value.values()):
        raise SpatialContractError(f"{where} dimensions must be positive integers")
    return value


def validate_contract(contract: dict[str, Any], vocabulary: dict[str, Any] | None = None) -> str:
    """Validate one artifact contract and return its canonical identity.

    The returned digest changes when geometry, coordinates, measures, field
    semantics, time support, operator order, or provenance changes. It is the
    identity callers stamp next to their data; it does not replace the hashes
    of the data sources named in ``provenance``.
    """
    vocab = vocabulary or declaration()
    if not isinstance(contract, dict):
        raise SpatialContractError("contract must be a mapping")
    for section in vocab["required_contract_sections"]:
        _need(contract, section, "contract")
    if contract["contract_version"] != vocab["contract_version"]:
        raise SpatialContractError(
            f"contract_version {contract['contract_version']!r} is not "
            f"{vocab['contract_version']!r}")
    if not isinstance(contract["artifact_id"], str) or not contract["artifact_id"].strip():
        raise SpatialContractError("artifact_id must be a non-empty string")

    supports = contract["supports"]
    if not isinstance(supports, list) or not supports:
        raise SpatialContractError("supports must be a non-empty list")
    support_by_id: dict[str, dict] = {}
    kinds = vocab["support_kinds"]
    for index, support in enumerate(supports):
        where = f"supports[{index}]"
        if not isinstance(support, dict):
            raise SpatialContractError(f"{where} must be a mapping")
        sid = _need(support, "id", where)
        if not isinstance(sid, str) or not sid:
            raise SpatialContractError(f"{where}.id must be a non-empty string")
        if sid in support_by_id:
            raise SpatialContractError(f"duplicate support id {sid!r}")
        kind = _need(support, "kind", where)
        if kind not in kinds:
            raise SpatialContractError(f"{where}.kind {kind!r} is not declared")
        _positive_dimensions(_need(support, "dimensions", where), where)
        for identity in ("geometry_sha256", "coordinates_sha256"):
            value = _need(support, identity, where)
            if not isinstance(value, str) or not _SHA256.fullmatch(value):
                raise SpatialContractError(f"{where}.{identity} must be a lowercase sha256")
        measures = _need(support, "native_measures", where)
        if not isinstance(measures, list) or not measures:
            raise SpatialContractError(f"{where}.native_measures must be a non-empty list")
        measure_kinds = set()
        for measure in measures:
            if not isinstance(measure, dict) or measure.get("kind") not in {
                    "area", "volume", "count"}:
                raise SpatialContractError(f"{where}.native_measures has an undeclared kind")
            if measure["kind"] in measure_kinds:
                raise SpatialContractError(f"{where} repeats native measure {measure['kind']!r}")
            if not measure.get("variable") or not measure.get("units"):
                raise SpatialContractError(f"{where}.native_measures need variable and units")
            measure_kinds.add(measure["kind"])
        if kind == "ocean_cgenie_grid" and not {"area", "volume"} <= measure_kinds:
            raise SpatialContractError(f"{where} ocean support requires area and volume")
        if kinds[kind].get("role") == "demographic_sample" and \
                support.get("spatially_representative") is not False:
            raise SpatialContractError(
                f"{where} is a demographic sample and must explicitly say "
                "spatially_representative: false")
        support_by_id[sid] = support

    fields = contract["fields"]
    if not isinstance(fields, list) or not fields:
        raise SpatialContractError("fields must be a non-empty list")
    field_by_name: dict[str, dict] = {}
    for index, field in enumerate(fields):
        where = f"fields[{index}]"
        if not isinstance(field, dict):
            raise SpatialContractError(f"{where} must be a mapping")
        name = _need(field, "name", where)
        if not isinstance(name, str) or not name or name in field_by_name:
            raise SpatialContractError(f"{where}.name must be non-empty and unique")
        semantics = _need(field, "semantics", where)
        if semantics not in vocab["field_semantics"]:
            raise SpatialContractError(f"{where}.semantics {semantics!r} is not declared")
        sid = _need(field, "support_id", where)
        if sid not in support_by_id:
            raise SpatialContractError(f"{where}.support_id {sid!r} is not in supports")
        spatial_dimensions = _need(field, "spatial_dimensions", where)
        if spatial_dimensions != list(support_by_id[sid]["dimensions"]):
            raise SpatialContractError(
                f"{where}.spatial_dimensions do not match support {sid!r}; "
                "shape equality is checked through named support dimensions")
        if not field.get("units"):
            raise SpatialContractError(f"{where}.units is required, including '1'")
        effective = _need(field, "effective_measure", where)
        if not isinstance(effective, dict) or effective.get("mode") not in \
                vocab["effective_measure_modes"]:
            raise SpatialContractError(f"{where}.effective_measure mode is not declared")
        if effective.get("measure_kind") not in {"area", "volume", "count"}:
            raise SpatialContractError(f"{where}.effective_measure.measure_kind is not declared")
        if effective["mode"] == "native_times_fraction" and not effective.get("fraction_field"):
            raise SpatialContractError(f"{where} needs effective_measure.fraction_field")
        if effective["mode"] == "explicit" and not all(
                effective.get(k) for k in ("variable", "units")):
            raise SpatialContractError(f"{where} explicit measure needs variable and units")
        aggregation = _need(field, "aggregation", where)
        if not isinstance(aggregation, dict):
            raise SpatialContractError(f"{where}.aggregation must be a mapping")
        if aggregation.get("operator") not in vocab["aggregation_operators"]:
            raise SpatialContractError(f"{where}.aggregation.operator is not declared")
        if aggregation.get("nonlinear_order") not in vocab["nonlinear_orders"]:
            raise SpatialContractError(f"{where}.aggregation.nonlinear_order is not declared")
        if not aggregation.get("operator_version"):
            raise SpatialContractError(f"{where}.aggregation.operator_version is required")
        for end in ("source_support", "destination_support"):
            if aggregation.get(end) not in support_by_id:
                raise SpatialContractError(f"{where}.aggregation.{end} is not in supports")
        if semantics == "categorical_label" and aggregation["operator"] not in {
                "identity", "categorical_majority"}:
            raise SpatialContractError(f"{where} categorical labels use a categorical operator")
        if semantics == "categorical_fraction" and aggregation["operator"] not in {
                "identity", "categorical_histogram", "area_weighted_fraction"}:
            raise SpatialContractError(f"{where} categorical fractions use a fraction operator")
        if semantics == "vector_component" and aggregation["operator"] not in {
                "identity", "vector_rotation_and_component_remap"}:
            raise SpatialContractError(f"{where} vectors require declared rotation/remap semantics")
        # The distribution pair is BIDIRECTIONAL, unlike the categorical pair
        # above. A quantile table is not a value of anything else, so the
        # operator implies the semantics as strongly as the semantics implies
        # the operator, and a field that declares one and not the other is a
        # cell mean wearing a distribution's name or the reverse.
        if semantics == "distribution_quantiles" and aggregation["operator"] not in {
                "identity", "area_weighted_distribution"}:
            raise SpatialContractError(
                f"{where} is a distribution and must carry the distribution operator; "
                "a mean is not a reduction of it, it is a different quantity")
        if aggregation["operator"] == "area_weighted_distribution" and \
                semantics != "distribution_quantiles":
            raise SpatialContractError(
                f"{where} carries the distribution operator, so its semantics is "
                "distribution_quantiles and not a single value per cell")
        if aggregation["operator"] == "area_weighted_distribution":
            _quantile_probabilities(field, where)
        if semantics == "vector_component" and field.get("vector_basis") not in \
                vocab["vector_bases"]:
            raise SpatialContractError(f"{where} vector_basis is not declared")
        field_by_name[name] = field

    for name, field in field_by_name.items():
        effective = field["effective_measure"]
        if effective["mode"] == "native_times_fraction":
            fraction = effective["fraction_field"]
            if fraction not in field_by_name:
                raise SpatialContractError(f"field {name!r} names missing fraction {fraction!r}")
            if field_by_name[fraction]["semantics"] not in {
                    "dimensionless_fraction", "categorical_fraction"}:
                raise SpatialContractError(f"field {name!r} effective fraction is not a fraction")
            if field_by_name[fraction]["support_id"] != field["support_id"]:
                raise SpatialContractError(f"field {name!r} fraction lives on another support")

    time = contract["time_support"]
    if not isinstance(time, dict) or time.get("semantics") not in vocab["time_semantics"]:
        raise SpatialContractError("time_support.semantics is not declared")
    if time["semantics"] != "static":
        for key in ("calendar", "absolute_time_units", "bounds_variable"):
            if not time.get(key):
                raise SpatialContractError(f"non-static time_support requires {key}")
        if time["semantics"].startswith("interval_") and time.get("interval_closure") not in {
                "[start,end)", "(start,end]"}:
            raise SpatialContractError("interval time_support needs an interval_closure")

    provenance = contract["provenance"]
    if not isinstance(provenance, dict) or not isinstance(
            provenance.get("source_artifacts"), list) or not provenance["source_artifacts"]:
        raise SpatialContractError("provenance.source_artifacts must be a non-empty list")
    if not isinstance(provenance.get("operator_versions"), dict) or not \
            provenance["operator_versions"]:
        raise SpatialContractError("provenance.operator_versions must be a non-empty mapping")
    return canonical_digest(contract)


def validate_conversion_assessment(
        assessment: dict[str, Any], source_contract: dict[str, Any],
        destination_contract: dict[str, Any],
        vocabulary: dict[str, Any] | None = None) -> str:
    """Validate one SPAT-10 conversion assessment and return its digest."""
    spec = vocabulary or conversion_declaration()
    if not isinstance(assessment, dict):
        raise SpatialContractError("conversion assessment must be a mapping")
    if assessment.get("assessment_version") != spec["assessment_version"]:
        raise SpatialContractError("conversion assessment version is not accepted")
    src_identity = validate_contract(source_contract)
    dst_identity = validate_contract(destination_contract)
    if assessment.get("source_contract_identity") != src_identity:
        raise SpatialContractError("source contract identity does not match its contract")
    if assessment.get("destination_contract_identity") != dst_identity:
        raise SpatialContractError("destination contract identity does not match its contract")

    def selected(contract: dict, sid: Any, end: str) -> dict:
        matches = [s for s in contract["supports"] if s["id"] == sid]
        if len(matches) != 1:
            raise SpatialContractError(f"{end}_support_id is not unique in its contract")
        return matches[0]

    src = selected(source_contract, assessment.get("source_support_id"), "source")
    dst = selected(destination_contract, assessment.get("destination_support_id"), "destination")
    if assessment.get("source_shape") != list(src["dimensions"].values()):
        raise SpatialContractError("source shape does not match the named support")
    if assessment.get("destination_shape") != list(dst["dimensions"].values()):
        raise SpatialContractError("destination shape does not match the named support")
    if assessment.get("source_coordinates_sha256") != src["coordinates_sha256"]:
        raise SpatialContractError("source coordinates do not match the named support")
    if assessment.get("destination_coordinates_sha256") != dst["coordinates_sha256"]:
        raise SpatialContractError("destination coordinates do not match the named support")
    if not assessment.get("operator_version"):
        raise SpatialContractError("operator_version is required")

    closures = assessment.get("closure")
    if not isinstance(closures, dict) or set(closures) != set(spec["closure_ledgers"]):
        raise SpatialContractError("closure ledger is incomplete or contains undeclared quantities")
    for name, row in closures.items():
        if not isinstance(row, dict) or row.get("status") not in spec["test_status"]:
            raise SpatialContractError(f"closure.{name} has no declared status")
        if row["status"] == "not_applicable":
            if not row.get("reason"):
                raise SpatialContractError(f"closure.{name} needs an inapplicability reason")
            continue
        if row.get("comparator") not in spec["comparators"]:
            raise SpatialContractError(f"closure.{name} comparator is not declared")
        if not all(isinstance(row.get(k), (int, float)) and not isinstance(row.get(k), bool)
                   for k in ("residual", "tolerance")):
            raise SpatialContractError(f"closure.{name} needs numeric residual and tolerance")
        passed = (row["residual"] <= row["tolerance"] if
                  row["comparator"] == "less_than_or_equal" else
                  row["residual"] >= row["tolerance"])
        if not passed:
            raise SpatialContractError(f"closure.{name} fails its declared tolerance")

    tests = assessment.get("tests")
    if not isinstance(tests, dict) or set(tests) != set(spec["required_tests"]):
        raise SpatialContractError("required conversion tests are incomplete")
    for name, row in tests.items():
        if not isinstance(row, dict) or row.get("status") not in spec["test_status"]:
            raise SpatialContractError(f"tests.{name} has no declared status")
        if row["status"] == "not_applicable":
            if not row.get("reason"):
                raise SpatialContractError(f"tests.{name} needs an inapplicability reason")
            continue
        if row.get("comparator") not in spec["comparators"]:
            raise SpatialContractError(f"tests.{name} comparator is not declared")
        if not all(isinstance(row.get(k), (int, float)) and not isinstance(row.get(k), bool)
                   for k in ("residual", "tolerance")):
            raise SpatialContractError(f"tests.{name} needs numeric residual and tolerance")
        passed = (row["residual"] <= row["tolerance"] if
                  row["comparator"] == "less_than_or_equal" else
                  row["residual"] >= row["tolerance"])
        if not passed:
            raise SpatialContractError(f"tests.{name} fails its declared tolerance")

    inventories = assessment.get("inventory")
    if not isinstance(inventories, dict) or set(inventories) != \
            set(spec["required_inventories"]):
        raise SpatialContractError("conversion inventory is incomplete")
    if any(not isinstance(value, list) for value in inventories.values()):
        raise SpatialContractError("every conversion inventory is a list, including empty ones")
    provenance = assessment.get("provenance")
    if not isinstance(provenance, dict) or not all(
            provenance.get(key) for key in
            ("source_artifact_sha256", "destination_geometry_sha256",
             "operator_sha256")):
        raise SpatialContractError("conversion provenance is incomplete")
    for key in ("source_artifact_sha256", "destination_geometry_sha256",
                "operator_sha256"):
        if not isinstance(provenance[key], str) or not _SHA256.fullmatch(provenance[key]):
            raise SpatialContractError(f"conversion provenance {key} is not a sha256")
    return canonical_digest(assessment)
