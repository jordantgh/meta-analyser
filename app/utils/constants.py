from enum import Enum, auto

class PageIdentity(Enum):
    SEARCH = auto()
    PARSED = auto()
    PRUNED = auto()

class Mode(Enum):
    BROWSING = auto()
    SEARCHING = auto()
    PROCESSING = auto()
    PRUNING = auto()