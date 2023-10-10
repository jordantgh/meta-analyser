from pyparsing import infixNotation, opAssoc, Keyword, Word, alphanums, QuotedString
import argparse
import os
import shutil

import re

class QueryTerm:
    def __init__(self, term):
        self.term = re.compile(re.escape(term[0]), re.IGNORECASE)

    def match(self, doc):
        return bool(self.term.search(doc))

class QueryAnd:
    def __init__(self, tokens):
        self.tokens = tokens[0]

    def match(self, doc):
        for token in self.tokens[::2]:
            if not token.match(doc):
                return False
        return True

class QueryOr:
    def __init__(self, tokens):
        self.tokens = tokens[0]

    def match(self, doc):
        for token in self.tokens[::2]:
            if token.match(doc):
                return True
        return False

class QueryNot:
    def __init__(self, tokens):
        self.query = tokens[0][1]

    def match(self, doc):
        return not self.query.match(doc)

def create_query_parser():
    AND = Keyword("AND")
    OR = Keyword("OR")
    NOT = Keyword("NOT")

    special_chars = "-_."
    alphanumeric_and_special = alphanums + special_chars

    term = (Word(alphanumeric_and_special) | QuotedString("'") | QuotedString('"')).setParseAction(QueryTerm)
    
    query_expr = infixNotation(term,
        [
            (NOT, 1, opAssoc.RIGHT, QueryNot),
            (AND, 2, opAssoc.LEFT, QueryAnd),
            (OR, 2, opAssoc.LEFT, QueryOr),
        ]
    )

    return query_expr

query_parser = create_query_parser()

def search(query, docs):
    query_expr = query_parser.parseString(query, parseAll=True)[0]
    return [doc[0] for doc in docs if query_expr.match(doc[1])] 

def load_documents(directory):
    docs = []
    for filename in os.listdir(directory):
        with open(os.path.join(directory, filename), 'r') as file:
            docs.append((filename, file.read()))
    return docs

def main():
    parser = argparse.ArgumentParser(description='Search plaintext files with a Boolean query.')
    parser.add_argument('-q', '--query', type=str, nargs='?', default=None, help='Query string')
    parser.add_argument('-qf', '--queryfile', type=str, nargs='?', default=None, help='File containing the query')
    parser.add_argument('targetdir', type=str, help='Directory to search')
    parser.add_argument('destdir', type=str, help='Destination directory')  # Add a new argument for the destination directory

    args = parser.parse_args()
    
    os.makedirs(args.destdir, exist_ok=True)

    if args.query is not None:
        query = args.query
    elif args.queryfile is not None:
        with open(args.queryfile, 'r') as file:
            query = file.read()
    else:
        print('Error: Either a query (-q) or a query file (-qf) must be provided.')
        return

    docs = load_documents(args.targetdir)
    matches = search(query, docs)

    for match in matches:
        print(match)
        # Copy each matched file to the destination directory
        shutil.copy(os.path.join(args.targetdir, match), os.path.join(args.destdir, match))


if __name__ == "__main__":
    main()
