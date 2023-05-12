# The primary function of this module is to convert a Duo fitting output to a
# new Duo input file. However, it has deliberately been written to be easily
# extensible, such that new actions to perform on an output file can easily be
# added.

import re

### Internal Duo object IDs ###

# All Duo objects are mapped to an internal ID number. New object types are 
# being added to Duo all the time to add a new object to the dictionary below 
# include its name string and its Duo internal ID number. You may need to update
# the ID numbers of other objects at the same time if their ID has changed.

object_type_ids = {
        "POTEN": 1,
        "POTENTIAL": 1,
        "SPIN-ORBIT": 2,
        "SPINORBIT": 2,
        "SPIN-ORBIT-X": 2,
        "L2": 3,
        "L**2": 3,
        "L^2": 3,
        "LXLY": 4,
        "LYLX": 4,
        "L+": 4,
        "L_+": 4,
        "LX": 4,
        "SPIN-SPIN": 5,
        "SPIN-SPIN-O": 6,
        "BOB-ROT": 7,
        "BOBROT" : 7,
        "SPIN-ROT": 8,
        "SPIN-ROTATION": 8,
        "DIABAT": 9,
        "DIABATIC": 9,
        "LAMBDA-OPQ": 10,
        "LAMBDA-P2Q": 11,
        "LAMBDAP2Q": 11,
        "LAMBDA-Q": 12,
        "LAMBDAQ": 12,
        "NAC": 13,
        "HFCC-BF": 21,
        "HFCC-A": 22,
        "HFCC-C": 23,
        "HFCC-D": 24,
        "HFCC-CI": 25,
        "HFCC-EQQ0": 26,
        "HFCC-EQQ2": 27,
        "QUADRUPOLE": 29,
        "ABINITIO": 30,
        "DIPOLE": 31,
        "TM": 32,
        "DIPOLE-MOMENT": 33,
        "DIPOLE-X": 34
        }


### File block classes ###

# Duo input files are composed of blocks that define objects or program arguments.
# All the information required to build a new input file can be obtained by 
# copying and replacing blocks of text from different parts of the output file.

# Class to store the position of a block of text within the file
class fileBlock:
    """Stores information about a block of text within a file.

    Attributes
    ----------
    lines : list[str]
        a list of strings corresponding to each line of text in the block
    start_line : int
        the line number of the start of the block in the file
    line_nums : list[int]
        a list of numbers corresponding to the position of each line of the block in the file
    """

    def __init__(self, start_line, init_lines=[]):
        """Create a fileBlock object

        Arguments
        ----------
        start_line : int
            the line number of the start of the block in the file
        init_lines : list[str]
            a list of strings corresponding to each line of text in the block

        """
        self.lines = init_lines
        self.start_line = start_line

    @property
    def line_nums(self): # return range of line numbers contained in block
        return list(range(self.start_line, self.start_line+len(self.lines)+1))

class objBlock(fileBlock):
    """Stores information about a block of text within a file representing a Duo object.

    Attributes
    ----------
    glob_param_line_nums : list[int]
        the global position in the file of the lines corresponding to the object's numerical parameters
    loc_param_line_nums : list[int]
        the local position relative to the start of the object block of the object's numerical parameters
    duo_type : str
        the Duo representation of the object (e.g EMO, MORSE, GRID, etc.)
    param_lines : list[str]
        the lines specifying the value of each numerical parameter
    """

    @property
    def glob_param_line_nums(self):
        for num, line in enumerate(self.lines):
            if line.split()[0].upper() == 'VALUES':
                glob_param_line_nums = range(self.start_line + num + 1, self.start_line + len(self.lines)-1)
                return list(glob_param_line_nums)
        return None
    
    @property
    def loc_param_line_nums(self):
        for num, line in enumerate(self.lines):
            if line.split()[0].upper() == 'VALUES':
                loc_param_line_nums = range(num + 1, len(self.lines))
                return list(loc_param_line_nums)
        return None

    @property
    def duo_type(self):
        for num, line in enumerate(self.lines):
            if line.split()[0].upper() == 'TYPE':
                return line.split()[1].upper()
        return None

    @property
    def param_lines(self):
        return [self.lines[i] for i in self.loc_param_line_nums]


### The actual output reader ###

class OutputReader:
    """Base class for processing Duo output file line by line.
    
    This base class performs no actions of its own. In order to modify an output
    file or generate a new input file a child class should be defined, that
    provides methods for handling the output.

    To add functionality:
    1. Create a new output reader object that inherits from the OutputReader class
    2. This new class must define a dictionary of triggers and methods
    3. The trigger should be a string to search for in the Duo output file.
    4. The method should be a function to pass the subsequent lines to when 
        the trigger is found
    5. The method must accept two arguments: the line number (int) and the line
        contents (str)
    5. Ensure these methods return 'False' until a stopping condition is met 
        e.g a stopping trigger line is found
    
    See IterationReader for an example.
    """

    def __init__(self):
        self.count_id  = 0
        self.waiting   = {}
        self.working   = {"triggers" : (self._check_triggers, {})} #starting state for the input reader
        self.triggers  = {}

    # Outline of file processing flow
    # -------------------------------
    # 1. Iterates over lines in the file and passes the contents of each line to any processing function currently 
    # in the working state (i.e those contained in `self.working'). If the function returns `True' then it is 
    # removed from the working list. 
    # 2. If the line matches one of triggers contained in `self.triggers' then the corresponding method is 
    # added to the waiting list (`self.waiting'). The function for checking each line for triggers is included
    # in the working list by default (`self._check_triggers`).
    # 3. After all working actions are processed, items in the waiting list are added to the working list 
    # before proceeding to the next line.

    def process_output(self, read_file):
        """Process the Duo output file line-by-line.
        
        Arguments
        ---------
        read_file : str
            The path to the Duo output file to be processed
        """
        # Initialise generator with the Duo output file
        with open(read_file, "r") as f: # read Duo output file line by line
            for num, line_ in enumerate(f):
                pop_ids = [] # empty list of actions to pop
                line = line_.rstrip('\n')
                for id_ in self.working: # perform waiting actions
                    f = self.working[id_][0]
                    kwargs = self.working[id_][1]
                    done = f(num, line, **kwargs) 
                    if done: # pop action if complete
                        pop_ids.append(id_)
                [self.working.pop(id_) for id_ in pop_ids]
                for id_ in self.waiting: # add waiting itemss
                    self.working[id_] = self.waiting[id_]
                self.waiting = {}

    def _check_triggers(self, num, line):
        """Checks the file line to see if it matches one of the stored trigger lines"""
        if line in self.triggers:
            for f in self.triggers[line]:
                self._new_waiting(f, {})
            return False
        else:
            return False
    
    def _new_waiting(self, func, kwargs):   
        self.count_id += 1
        self.waiting[str(self.count_id)] = (func, kwargs)

class IterationReader(OutputReader):
    """Duo output processor for updating Duo object parameters based on the fitted values from a given iteration."""

    def __init__(self):
        super().__init__()
        self.input_transcript = fileBlock(39)
        self.input_objects    = {}
        self.iter_objects     = []
        self.triggers = {
            "(Transcript of the input --->)" : [self._read_input_transcript],
            "Parameters:" : [self._add_new_iteration, self._read_fitting_parameters]
        } # lines in Duo output and the action(s) they trigger

    def genfromit(self, it_num=-1):
        """Generate a new fitting input from the fitted parameters from a given iteration number.
        
        Arguments
        ---------
        it_num : int
            The iteration number to take the fitted parameters from (default is last iteration)
        """
        for num, line in enumerate(self.input_transcript.lines): #read Duo output file line by line
            for obj_id, obj in self.iter_objects[it_num].items():
                if (num not in self.input_objects[obj_id].glob_param_line_nums):
                    continue
                else:
                    inpline = line.split()
                    outline = obj.param_lines[num - self.input_objects[obj_id].glob_param_line_nums[0]].split()
                    
                    preval = inpline[1]
                    if self.input_objects[obj_id].duo_type == "GRID":
                        pass #grid type potentials do not change during fitting
                    elif len(inpline) > 2 and inpline[2] == "fit":
                        if len(inpline) > 3 and inpline[3] == "link":
                            line = "{0:11} {1: .14E} fit link {2} {3} {4} ({5: .14E})".format(outline[0], float(outline[1]), *outline[3:6], float(preval))
                        else:
                            line = "{0:11} {1: .14E} fit            ({2: .14E})".format(outline[0], float(outline[1]), float(preval))
                    else:
                        line = "{0:11} {1: .14E}                ({2: .14E})".format(outline[0], float(outline[1]), float(preval))
            yield line
    
    def _read_input_transcript(self, num, line):
        """Detects the presence of a Duo object in the input section, initialises an `objBlock' and adds the associated method to the waiting list"""
        if line == "(<--- End of the input.)":
            return True
        else:
            self.input_transcript.lines.append(line)
            if any(re.search(rf"^\s*{obj}", line.upper()) for obj in object_type_ids):
                obj_id = self._determine_id(line)
                self.input_objects[obj_id] = objBlock(num-self.input_transcript.start_line, init_lines=[line])
                self._new_waiting(self._read_input_objects, {'obj_id' : obj_id})
            return False

    def _read_fitting_parameters(self, num, line):
        """Detects the presence of a Duo object in the fitting section, initialises an `objBlock' and adds the associated method to the waiting list"""
        if line == "Fitted paramters (rounded):":
            return True
        else:
            if any(re.search(rf"^\s*{obj}", line.upper()) for obj in object_type_ids):
                obj_id = self._determine_id(line)
                self.iter_objects[-1][obj_id] = objBlock(num, init_lines=[line])
                self._new_waiting(self._read_fitting_objects, {'obj_id' : obj_id})
    
    def _add_new_iteration(self, num, line):
        """Adds a new object to store the results of an iteration in the Duo fitting procedure"""
        self.iter_objects.append({})
        return True
    
    def _read_input_objects(self, num, line, obj_id=None):
        """Appends the current line to the list of lines describing a Duo object in the original input until the end of the Duo object block is reached"""
        self.input_objects[obj_id].lines.append(line)
        if line.strip().upper() == "END":
            return True
        else:
            return False
    
    def _read_fitting_objects(self, num, line, obj_id=None):
        """Appends the current line to the list of lines describing a Duo object in the fitting output until the end of the Duo object block is reached"""
        self.iter_objects[-1][obj_id].lines.append(line)
        if line.strip().upper() == "END":
            return True
        else:
            return False

    def _determine_id(self, line):
        """Generate a unique ID tuple in order to track all the fitted items.

        Arguments
        ----------
        line : str
            A file line corresponding to the start of the Duo object

        Returns
        -------
        obj_id : tuple(Bool, int, int, int)
            A tuple that uniquely identifies the Duo object. Contains a Boolean flag that indicates whether
            the object type is proceeded by an ABINITIO flag (see Duo manual); an integer number enumerating
            the object type (e.g POTEN, DIPOLE, SPIN-ORBIT); and two IDs indicating the index of the
            POTEN objects coupled by a non-POTEN object, or the POTEN index itself.
        """
        keys = line.upper().split()
        obj_type = object_type_ids[keys[0]]
        # ABINITIO objects are special case of another Duo object
        if obj_type == object_type_ids["ABINITIO"]:
            abinitio = True
            try:
                obj_type = object_type_ids[keys[1]]
            except KeyError as error:
                print("Found unknown abinitio object. This usually occurs when a new object type has been added to Duo. Try adding the object and it's internal Duo ID to the object_type_ids dictionary.")
                raise error
            keys = keys[1:]    
        else:
            abinitio = False
        # POTEN objects have a single state ID, couplings have two state IDs
        if obj_type == object_type_ids["POTEN"]:
            l_id = r_id = int(keys[1])
        else:
            l_id, r_id = int(keys[1]), int(keys[2])
        obj_id = (abinitio, obj_type, l_id, r_id)
        return obj_id
    

### Command line argument parser ###

import argparse

# Define runtime arguments that can be passed to the script
parser = argparse.ArgumentParser(description="Duo fitting input iterator")
parser.add_argument(
    'input', metavar='duo_output.out', type=str,
    help="Reference file to generate new input from."
    )
parser.add_argument(
    '-f', '--file', metavar='my_input.inp', type=str,
    help="Name of the Duo '.inp' file to write output to, if not console." 
)
parser.add_argument(
    '-n', '--number', metavar='n', type=int, default=-1,
    help="The iteration number (zero indexed) from which to take the fitted parameters. Defaults to the last iteration. Can also be specified from the end using negative numbers, i.e -n -2 uses the penultimate iteration."
)
args = parser.parse_args() #parse arguments

if __name__ == "__main__":
    gen = IterationReader()
    gen.process_output(args.input)

    fout = None
    if args.file is not None:
        fout = open(args.file, 'w+')
    for line in gen.genfromit(it_num=args.number):
        print(line, file=fout)
    if args.file is not None:
        fout.close()

