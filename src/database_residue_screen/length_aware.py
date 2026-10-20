"""Length-aware score retention for the fixed blastn-short scoring contract.

The length adjustment is a Python translation of NCBI's public-domain
BLAST_ComputeLengthAdjustment (blast_stat.c, BLAST+ 2.17.0). Runtime XML
probes check the statistics, local search space and exact witness recovery.
No study database dimensions or panel lengths are built into this module.
"""
import math
import os
import random
import subprocess
import tempfile
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

from .contracts import ContractError
from .example import save_fasta


def length_adjustment(length, bases, sequences, kappa, entropy):
    if not (0 < length < 2**31 and 0 < bases < 2**63 and 0 < sequences < 2**31
            and kappa > 0 and entropy > 0):
        raise ContractError("Invalid BLAST statistical inputs")
    m, n, count = float(length), float(bases), float(sequences)
    mb = m * count + n
    c = n * m - max(m, n) / kappa
    if c < 0:
        return 0
    upper = 2 * c / (mb + math.sqrt(max(0, mb * mb - 4 * count * c)))
    lower, next_value, converged = 0., 0., False
    for iteration in range(1, 21):
        value = next_value
        space = (m - value) * (n - count * value)
        proposed = (math.log(kappa) + math.log(space)) / entropy
        if proposed >= value:
            lower = value
            if proposed - lower <= 1:
                converged = True
                break
            if lower == upper:
                break
        else:
            upper = value
        next_value = proposed if lower <= proposed <= upper else upper if iteration == 1 else (lower + upper) / 2
    adjustment = int(lower)
    if converged:
        value = math.ceil(lower)
        if value <= upper and (math.log(kappa) + math.log((m - value) * (n - count * value))) / entropy >= value:
            adjustment = value
    return adjustment


def minimum_tolerant_score(length):
    """All ungapped spans ceil(.90 L)..L, at least 95% identical; +1/-3."""
    if length < 1:
        raise ContractError("Query length must be positive")
    return min(span - 4 * (span // 20) for span in range((9 * length + 9) // 10, length + 1))


def chance_bound(length, bases, endpoint):
    """Uniform iid both-strand occurrence expectation (tolerant: union bound).

    Count each allowed query span and substitution pattern, including overlapping
    spans. This is a conservative model bound, not an empirical contamination
    probability; repeats/composition require separate matched-control assessment.
    Statistical -searchsp never changes the actual database opportunity count.
    """
    if endpoint == "strict":
        return 2 * bases * math.exp(-length * math.log(4))
    return 2 * bases * sum((length - span + 1) * math.exp(-span * math.log(4)) *
                          sum(math.comb(span, errors) * 3**errors for errors in range(span // 20 + 1))
                          for span in range((9 * length + 9) // 10, length + 1))


def testability(retainable, expected_chance, informative_limit=.05, testable_limit=.01):
    if not (0 < testable_limit <= informative_limit < math.inf):
        raise ContractError("Require 0 < testable chance limit <= informativeness limit < infinity")
    if not retainable or expected_chance > informative_limit:
        return "censored_untestable"
    return "testable" if expected_chance <= testable_limit else "chance_limited"


def testability_reason(score_retainable, witness_recovered, expected_chance, informative_limit=.05, testable_limit=.01):
    """Explain a class; `censored_untestable` covers score censoring AND chance saturation.

    A retainable 19-nt exact query is untestable only because chance matches are
    expected, not because BLAST censors it; reports must not conflate the two.
    """
    state = testability(score_retainable and witness_recovered, expected_chance, informative_limit, testable_limit)
    if state == "testable":
        return "chance_bound_within_testable_limit"
    if state == "chance_limited":
        return "chance_bound_between_limits"
    return "+".join(name for name, failed in (("score_not_retainable", not score_retainable),
                                               ("no_exact_witness", not witness_recovered),
                                               ("chance_bound_above_informativeness_limit", expected_chance > informative_limit))
                    if failed)


def score_calibration(length, bases, sequences, statistics, evalue=10., mode="fixed", maximum_evalue=None, searchsp=None,
                      informative_limit=.05, testable_limit=.01):
    if mode not in {"adaptive", "fixed"} or not 0 < evalue < math.inf:
        raise ContractError("Invalid length-aware threshold policy")
    if maximum_evalue is not None and not 0 < maximum_evalue < math.inf:
        raise ContractError("Maximum E-value must be finite and positive")
    kappa, lam, entropy = (statistics[name] for name in ("kappa", "lambda_value", "entropy"))
    adjustment = length_adjustment(length, bases, sequences, kappa, entropy)
    space = searchsp if searchsp is not None else max(1, bases - sequences * adjustment) * (length - adjustment)
    if space < 1:
        raise ContractError("Effective search space must be positive")
    exact = kappa * space * math.exp(-lam * length)
    minimum_score = minimum_tolerant_score(length)
    tolerant = kappa * space * math.exp(-lam * minimum_score)
    threshold = max(evalue, 1.1 * tolerant) if mode == "adaptive" else evalue
    if maximum_evalue is not None:
        threshold = min(threshold, maximum_evalue)
    # Below 4 nt BLAST cannot seed. Longer endpoints always contain a 7-mer
    # under 95% identity; short queries use a smaller legal word.
    word = max(4, min(7, length))
    strict_retained = length >= word and exact <= threshold
    tolerant_retained = length >= word and tolerant <= threshold
    strict_chance, tolerant_chance = (chance_bound(length, bases, endpoint) for endpoint in ("strict", "tolerant"))
    return dict(query_length=length, database_bases=bases, database_sequences=sequences,
                length_adjustment=adjustment, effective_search_space=space, **statistics,
                exact_score=length, exact_evalue=exact, minimum_tolerant_score=minimum_score,
                tolerant_evalue=tolerant, evalue=threshold, word_size=word,
                strict_testable=strict_retained, tolerant_fully_testable=tolerant_retained,
                strict_score_retainable=strict_retained, tolerant_score_retainable=tolerant_retained,
                partially_score_censored=strict_retained and not tolerant_retained,
                strict_expected_chance_bound=strict_chance, tolerant_expected_chance_bound=tolerant_chance,
                strict_testability=testability(strict_retained, strict_chance, informative_limit, testable_limit),
                tolerant_testability=testability(tolerant_retained, tolerant_chance, informative_limit, testable_limit),
                # Runtime witnesses are added by calibrate_lengths, which recomputes both reasons.
                strict_testability_reason=testability_reason(strict_retained, True, strict_chance, informative_limit, testable_limit),
                tolerant_testability_reason=testability_reason(tolerant_retained, True, tolerant_chance, informative_limit, testable_limit),
                informativeness_limit=informative_limit, testable_chance_limit=testable_limit,
                diagnostic_chance_saturated_output=mode == "adaptive",
                # KA E counts chance HSPs scoring at least S, not accepted constructs.
                expected_random_hsps_at_exact_score=exact,
                expected_random_hsps_at_tolerant_score=tolerant,
                expected_random_hsps_at_search_threshold=threshold,
                iid_uniform_exact_occurrences_both_strands=2 * bases * math.exp(-length * math.log(4)))


def engine_statistics(binary, lengths):
    """Calibrate against a tiny generated DB, per query, using the actual binary.

    stdin avoids persistent temporary query paths in provenance. Returned values
    contain no temporary paths or timestamps so resumes remain deterministic.
    """
    generator = random.Random(19019)
    requested_lengths = sorted(set(lengths))
    probes = {length: ''.join(generator.choices('ACGT', k=length)) for length in sorted(set(lengths) | {32}) if length >= 4}
    subjects = {"probe": 'N'.join('N' * 50 + seq + 'N' * 50 for seq in probes.values())}
    db_bases = len(subjects['probe'])
    with tempfile.TemporaryDirectory(prefix="residue-calibration-") as temporary:
        root = Path(temporary)
        save_fasta(root / "subject.fa", subjects)
        builder = str(Path(binary).with_name("makeblastdb"))
        subprocess.run([builder, "-in", str(root / "subject.fa"), "-dbtype", "nucl", "-out", str(root / "db")],
                       capture_output=True, text=True, check=True, timeout=30)
        checked, reference = [], None
        for length, sequence in probes.items():
            command = [binary, "-db", str(root / "db"), "-task", "blastn-short", "-word_size", str(max(4, min(7, length))),
                       "-evalue", "1e12", "-dust", "no", "-soft_masking", "false", "-reward", "1", "-penalty", "-3",
                       "-gapopen", "5", "-gapextend", "2", "-outfmt", "5", "-num_threads", "1"]
            result = subprocess.run(command, input=f">q\n{sequence}\n", capture_output=True, text=True, check=True,
                                    env=dict(os.environ, LC_ALL="C", BLAST_USAGE_REPORT="false"), timeout=30)
            root_xml = ET.fromstring(result.stdout)
            stat = root_xml.find('.//Statistics')
            if stat is None:
                raise ContractError("BLAST calibration did not emit XML statistics")
            observed = dict(kappa=float(stat.findtext('Statistics_kappa')), lambda_value=float(stat.findtext('Statistics_lambda')),
                            entropy=float(stat.findtext('Statistics_entropy')))
            if reference is not None and any(not math.isclose(observed[k], reference[k], rel_tol=1e-10) for k in observed):
                raise ContractError("Length-only statistical calibration is invalid for this engine")
            reference = observed
            row = score_calibration(length, db_bases, 1, observed)
            actual_space = int(stat.findtext('Statistics_eff-space'))
            if row['effective_search_space'] != actual_space:
                raise ContractError("Length adjustment disagrees with the installed BLAST engine")
            for hsp in root_xml.findall('.//Hsp'):
                expected = observed['kappa'] * actual_space * math.exp(-observed['lambda_value'] * float(hsp.findtext('Hsp_score')))
                if not math.isclose(expected, float(hsp.findtext('Hsp_evalue')), rel_tol=1e-5, abs_tol=1e-300):
                    raise ContractError("Score/E-value calibration disagrees with BLAST")
            recovered = any(int(h.findtext('Hsp_score')) == length and int(h.findtext('Hsp_align-len')) == length
                            for h in root_xml.findall('.//Hsp'))
            checked.append(dict(length=length, exact_witness_recovered=recovered, local_effective_search_space=actual_space))
    by_length = {row['length']: row for row in checked}
    return reference, [by_length.get(length, dict(length=length, exact_witness_recovered=False, local_effective_search_space=0)) for length in requested_lengths]


def score_retained(row, endpoint):
    """Score-only retention of a calibration row, before runtime witness checks."""
    evalue = row['exact_evalue'] if endpoint == 'strict' else row['tolerant_evalue']
    return row['query_length'] >= row['word_size'] and evalue <= row['evalue']


def verify_calibration_row(row):
    """Recompute every derived field that can be recomputed from the row's own inputs.

    Guards hand-edited or migrated calibrations: classes, reasons, chance bounds and
    score E-values must follow from length, database size, statistics and threshold.
    """
    def close(left, right):
        return math.isclose(left, right, rel_tol=1e-9, abs_tol=1e-300)
    length, limits = row['query_length'], (row['informativeness_limit'], row['testable_chance_limit'])
    if 'minimum_tolerant_score' in row and row['minimum_tolerant_score'] != minimum_tolerant_score(length):
        raise ContractError("Calibration minimum tolerant score is inconsistent")
    if {'kappa', 'lambda_value', 'effective_search_space'} <= row.keys():
        for field, score in (('exact_evalue', length), ('tolerant_evalue', minimum_tolerant_score(length))):
            if not close(row[field], row['kappa'] * row['effective_search_space'] * math.exp(-row['lambda_value'] * score)):
                raise ContractError("Calibration E-value is inconsistent with its statistics")
    for endpoint, retain in (("strict", "strict_testable"), ("tolerant", "tolerant_fully_testable")):
        chance = row[endpoint + '_expected_chance_bound']
        if 'database_bases' in row and not close(chance, chance_bound(length, row['database_bases'], endpoint)):
            raise ContractError("Calibration chance bound is inconsistent with the database size")
        if row[retain] != (score_retained(row, endpoint) and row['exact_witness_recovered']):
            raise ContractError("Calibration retention flag is inconsistent with its scores and witness")
        if row[endpoint + '_testability'] != testability(row[retain], chance, *limits):
            raise ContractError("Testability class disagrees with declared calibration")
        reason = endpoint + '_testability_reason'
        if reason in row and row[reason] != testability_reason(score_retained(row, endpoint), row['exact_witness_recovered'], chance, *limits):
            raise ContractError("Testability reason disagrees with declared calibration")
    return row


def calibrate_lengths(binary, lengths, metadata, settings, mode="fixed", maximum_evalue=None,
                      informative_limit=.05, testable_limit=.01):
    statistics, witnesses = engine_statistics(binary, lengths)
    rows = [score_calibration(length, metadata['base_count'], metadata['sequence_count'], statistics,
                              settings['evalue'], mode, maximum_evalue, settings.get('searchsp'),
                              informative_limit, testable_limit) for length in sorted(set(lengths))]
    for row, witness in zip(rows, witnesses):
        row['exact_witness_recovered'] = witness['exact_witness_recovered']
        row['strict_testable'] &= witness['exact_witness_recovered']
        row['tolerant_fully_testable'] &= witness['exact_witness_recovered']
        for endpoint, key in (("strict", "strict_testable"), ("tolerant", "tolerant_fully_testable")):
            row[endpoint + '_score_retainable'] = row[key]
            row[endpoint + '_testability'] = testability(row[key], row[endpoint + '_expected_chance_bound'], informative_limit, testable_limit)
            row[endpoint + '_testability_reason'] = testability_reason(
                score_retained(row, endpoint), row['exact_witness_recovered'], row[endpoint + '_expected_chance_bound'],
                informative_limit, testable_limit)
        if any(row[e + '_testability'] != 'testable' for e in ('strict', 'tolerant')):
            warnings.warn(f"Length {row['query_length']}: {row['strict_testability']}/{row['tolerant_testability']} at E={row['evalue']:g}", RuntimeWarning)
    if mode == "adaptive":
        warnings.warn("DIAGNOSTIC ONLY: raised thresholds produce chance-saturated output; testability exclusions remain", RuntimeWarning)
    return dict(method="NCBI_length_adjustment_runtime_XML_checked_testability_v2", threshold_mode=mode,
                informativeness_limit=informative_limit, testable_chance_limit=testable_limit,
                chance_model="uniform_iid_both_strands_tolerant_union_bound_not_empirical",
                diagnostic_chance_saturated_output=mode == "adaptive",
                maximum_evalue=maximum_evalue, runtime_witnesses=witnesses, lengths=rows,
                limitation="Score retainability and isolated exact witnesses; heuristic HSP representation and target caps can still lose hits")
