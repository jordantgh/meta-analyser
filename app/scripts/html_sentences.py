from difflib import SequenceMatcher
from nltk.tokenize import sent_tokenize
from unicodedata import normalize


class InsertMatcher(SequenceMatcher):
    def get_opcodes(self) -> 'list[tuple]':
        opcodes = super().get_opcodes()
        new_opcodes = []
        for tag, i1, i2, j1, j2 in opcodes:
            if tag == 'replace':
                new_opcodes.append(('insert', i2, i2, j1, j2))
            else:
                new_opcodes.append((tag, i1, i2, j1, j2))
        return new_opcodes


def get_index_mappings(sentences, text) -> 'list[tuple[int, int]]':
    mappings = []
    last_idx = 0
    for sentence in sentences:
        # search from the last index to handle duplicate sentences
        start_idx = text.find(sentence, last_idx)
        if start_idx == -1:
            raise ValueError(
                "The sentence doesn't seem to match with the original text.")
        end_idx = start_idx + len(sentence)
        last_idx = end_idx
        mappings.append((start_idx, end_idx))

    return mappings


def filter_opcodes(
        start_idx: 'int',
        end_idx: 'int',
        opcodes: 'list[tuple]'
) -> 'list[tuple]':
    filtered_opcodes = []
    for tag, i1, i2, j1, j2 in opcodes:
        if i2 <= start_idx or i1 >= end_idx:
            continue
        filtered_opcodes.append(
            (tag, max(i1, start_idx) - start_idx, min(i2, end_idx) - start_idx, j1, j2))
    return filtered_opcodes


def apply_opcodes(
    sentence: 'str',
    filtered_opcodes: 'list[tuple]',
    original_html: 'str'
) -> 'str':
    new_text = []
    i1 = 0
    for tag, i1, i2, j1, j2 in filtered_opcodes:
        new_text.append(sentence[i1:i2])
        if tag == 'insert':
            new_text.append(original_html[j1:j2])

    new_sentence = ''.join(new_text)
    return new_sentence


def get_html_sentence(
    sentence: 'str',
    start_idx: 'int',
    end_idx: 'int',
    opcodes: 'list[tuple]',
    original_html: 'str'
) -> 'str':
    filtered_opcodes = filter_opcodes(start_idx, end_idx, opcodes)
    new_sentence = apply_opcodes(sentence, filtered_opcodes, original_html)

    return new_sentence


def get_sentences(text: 'str') -> 'list[str]':
    normalized_text = normalize('NFC', text)  # Normalize the text
    sentences = sent_tokenize(normalized_text)
    sentence_mappings = get_index_mappings(sentences, normalized_text)
    return sentences, sentence_mappings
