from abc import ABC, abstractmethod

class Model(ABC):
    """Base template of a Model class. Must implement a solve method!"""
    
    @abstractmethod
    def solve(self, *args, **kwargs):
        pass
    