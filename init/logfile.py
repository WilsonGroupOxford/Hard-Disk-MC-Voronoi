"""Simple log file class"""
from datetime import datetime
import sys

class Logfile:


    def __init__(self,name='logfile',**kwargs):
        """Open log file and write header"""

        # Get formatting options
        indent = kwargs.get('indent_len',4)
        dash = kwargs.get('dash_len',60)
        self.indent = ' '*indent
        self.dash = '-'*dash + '\n'

        # Print header
        self.file = open('./{}.log'.format(name),'w')
        self.time_start = datetime.now()
        self.file.write('Process begun: {}  \n'.format(self.time_start.strftime('%H:%M:%S')))
        self.dashed_line()


    def close(self,termination='normal'):
        """Close log file"""

        # Print footer
        self.time_end = datetime.now()
        time_delta = (self.time_end - self.time_start).total_seconds()
        self.file.write('Process completed with {} termination: {}  \n'.format(termination,self.time_end.strftime('%H:%M:%S')))
        self.file.write('Total run time: {} mins {} secs  \n'.format(*divmod(time_delta,60)))
        self.file.close()


    def __call__(self,line,indent=0,dash=False):
        """Write function as call for convenience"""

        line = '{}{}  \n'.format(indent*self.indent,line)
        self.file.write(line)
        if dash:
            self.dashed_line()


    def write(self,line,indent=0,dash=False):
        """Write line"""

        line = '{}{}  \n'.format(indent*self.indent,line)
        self.file.write(line)
        if dash:
            self.dashed_line()


    def error(self,message,type='critical'):
        """Write error message and kill process if neccessary"""

        if type=='critical':
            self.file.write('Critical error: {}  \n'.format(message))
            self.dashed_line()
            self.close(termination='error')
            sys.exit()


    def dashed_line(self):
        """Write dashed line"""

        self.file.write(self.dash)


if __name__ == '__main__':

    test = Logfile(name='test')
    test.write('Test Header')
    test.write('Test item 1',indent=1)
    test.write('Test item 2',indent=2)
    test.write('Test item 3',indent=2,dash=True)
    test.close()
