from itertools import product


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _witness(query_length, query_offset, subject_offset, span, mismatches, strand):
    return dict(qstart=query_offset + 1 if strand == "+" else query_length - query_offset,
                qend=query_offset + span if strand == "+" else query_length - query_offset - span + 1,
                sstart=subject_offset + 1, send=subject_offset + span, strand=strand,
                span=span, matches=span - mismatches, mismatches=mismatches)


def sequence_witnesses(query: str, subject: str) -> dict:
    if not query or set(query) - set("ACGT") or set(subject) - set("ACGTN"):
        raise ValueError("Fixture truth requires an ACGT query and ACGTN subject")
    length = len(query)
    minimum_span = (9 * length + 9) // 10
    orientations = (("+", query), ("-", reverse_complement(query)))
    for strand, oriented in orientations:
        position = subject.find(oriented)
        if position >= 0:
            exact = _witness(length, 0, position, length, 0, strand)
            return {"strict": dict(exact), "tolerant": dict(exact)}
    result = {"strict": None, "tolerant": None}
    if len(subject) < minimum_span:
        return result
    minimum_run = min((span - span // 20 + span // 20) // (span // 20 + 1)
                      for span in range(minimum_span, length + 1))
    for strand, oriented in orientations:
        if not any(oriented[offset:offset + minimum_run] in subject for offset in range(length - minimum_run + 1)):
            continue
        if length < 20:
            for span in range(length, minimum_span - 1, -1):
                for query_offset in range(length - span + 1):
                    position = subject.find(oriented[query_offset:query_offset + span])
                    if position >= 0:
                        result["tolerant"] = _witness(length, query_offset, position, span, 0, strand)
                        return result
            continue
        diagonals = []
        for diagonal in range(minimum_span - length, len(subject) - minimum_span + 1):
            query_offset, subject_offset = max(0, -diagonal), max(0, diagonal)
            overlap = min(length - query_offset, len(subject) - subject_offset)
            prefix = [0]
            for offset in range(overlap):
                prefix.append(prefix[-1] + (oriented[query_offset + offset] != subject[subject_offset + offset]))
            diagonals.append((query_offset, subject_offset, prefix))
        for span in range(length, minimum_span - 1, -1):
            for query_offset, subject_offset, prefix in diagonals:
                for offset in range(len(prefix) - span):
                    errors = prefix[offset + span] - prefix[offset]
                    if 20 * errors <= span:
                        result["tolerant"] = _witness(length, query_offset + offset, subject_offset + offset, span, errors, strand)
                        return result
    return result


def precompute_truth(queries: dict[str, str], subjects: dict[str, str]) -> tuple[dict, list]:
    truth = {"strict": set(), "tolerant": set()}
    witnesses = []
    for query_id, subject_id in product(sorted(queries), sorted(subjects)):
        for endpoint, witness in sequence_witnesses(queries[query_id], subjects[subject_id]).items():
            if witness is not None:
                truth[endpoint].add((query_id, subject_id))
                witnesses.append(dict(query_id=query_id, subject_id=subject_id, endpoint=endpoint, **witness))
    return {endpoint: frozenset(pairs) for endpoint, pairs in truth.items()}, witnesses


def independent_hsp_truth(row: dict, query: str, subject: str) -> dict:
    query_start, query_end = int(row["qstart"]), int(row["qend"])
    subject_start, subject_end = int(row["sstart"]), int(row["send"])
    query_low, query_high = sorted((query_start, query_end))
    subject_low, subject_high = sorted((subject_start, subject_end))
    if not (1 <= query_low <= query_high <= len(query) and 1 <= subject_low <= subject_high <= len(subject)):
        raise ValueError("Raw HSP coordinates leave the frozen sequence")
    aligned_query, aligned_subject = row["qseq"], row["sseq"]
    expected_query = query[query_low - 1:query_high]
    expected_subject = subject[subject_low - 1:subject_high]
    if query_start > query_end:
        expected_query = reverse_complement(expected_query)
    if subject_start > subject_end:
        expected_subject = reverse_complement(expected_subject)
    if aligned_query.replace("-", "") != expected_query or aligned_subject.replace("-", "") != expected_subject:
        raise ValueError("Raw aligned strings disagree with the frozen sequence")
    if len(aligned_query) != len(aligned_subject):
        raise ValueError("Raw alignment columns disagree")
    span = query_high - query_low + 1
    ungapped = "-" not in aligned_query + aligned_subject and len(aligned_query) == span == subject_high - subject_low + 1
    exact = sum(query_base == subject_base for query_base, subject_base in zip(aligned_query, aligned_subject))
    return dict(strict=ungapped and span == len(query) and exact == span,
                tolerant=ungapped and 10 * span >= 9 * len(query) and 20 * exact >= 19 * span)