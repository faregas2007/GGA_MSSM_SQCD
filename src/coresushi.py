# coresushi.py
# include sqcd c factor
# bidirectional control
from renschemes import *
from sigma import *
from include import *

# invoke: where all commands will be executed based on some condition.
class invoke:

    def __init__(self):
        self._commands = {}

    def register(self, command_name, command):
        """Register commands in Invoker """
        self._commands[command_name] = command
    
    def execute(self, command_name):
        """Execute any registered commands"""
        if(command_name in self._commands.keys()):
            self._commands[command_name].execute()
        else:
            print(f"Command[{command_name}] not recognised")


# Receiver --> in the imported modules.

# submit multiple jobs
class cluster(Command):
    def __init__(self, var: cluster):
        self._var = cluster()
    
    def execute(self):
        self.var.to_file()

# evaluate the cross-section
class distribution(Command):
    def __init__(self, var: sigma):
        self._var = sigma()
    
    def execute(self):
        self._var.to_file()

# doing the renormlize
class sushi_renormalize(Command):
    def __init__(self, var: renschemes)
        self._var = renschemes()

    def execute(self):
        self._var.to_file()

if __name__ == "__main__":
    ren = renschemes(inputs, einital, renormalizeSM, quarkhiggscoup, squarkhiggscoupCMSSM, renormalize)    

    command1 = sushi_renormalize(ren)
    
    invoker = invoker()
    invoker.register("ren", command1)

    invoker.execute("ren")