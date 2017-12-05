from optparse import OptionParser

from as_python.samseg.command_arguments import parse_args


class MockingParser(OptionParser):
    def __init__(self):
        super().__init__()
        self.error_message_received = None

    def get_prog_name(self):
        return "MockingParser"

    def error(self, msg):
        self.error_message_received = msg


class TestParseArgs:
    def setup(self):
        self.parser = None

    def parse_command_line(self, command_line):
        self.parser = MockingParser()
        args = command_line.split()
        return parse_args(args, parser=self.parser)

    def parse_minimum_command_line(self):
        return self.parse_command_line('-iread_here')

    def test_no_verbose_arg(self):
        recipe = self.parse_minimum_command_line()
        assert not recipe.verbose

    def test_short_verbose_arg(self):
        recipe = self.parse_command_line('-v')
        assert recipe.verbose

    def test_long_verbose_arg(self):
        recipe = self.parse_command_line('--verbose')
        assert recipe.verbose

    def test_no_exvivo_arg(self):
        recipe = self.parse_minimum_command_line()
        assert not recipe.exvivo

    def test_exvivo_arg(self):
        recipe = self.parse_command_line('--exvivo')
        assert recipe.exvivo

    def test_default_thread_count(self):
        recipe = self.parse_minimum_command_line()
        assert 1 == recipe.threads

    def test_thread_count(self):
        recipe = self.parse_command_line('--threads 5')
        assert 5 == recipe.threads

    def test_no_input(self):
        recipe = self.parse_command_line('')
        assert self.parser.error_message_received is not None

    def test_input_folder_short_flag(self):
        recipe = self.parse_command_line('-isome_folder')
        assert ['some_folder'] == recipe.image_file_names

    def test_input_folder_long_flag(self):
        recipe = self.parse_command_line('--input some_folder')
        assert ['some_folder'] == recipe.image_file_names

    def test_output_folder_short_flag(self):
        recipe = self.parse_command_line('-osome_folder')
        assert 'some_folder' == recipe.output

    def test_output_folder_long_flag(self):
        recipe = self.parse_command_line('--output some_folder')
        assert 'some_folder' == recipe.output

    def test_no_missing_structures(self):
        recipe = self.parse_minimum_command_line()
        assert [] == recipe.missing_structures

    def test_one_missing_structure_short(self):
        recipe = self.parse_command_line('-mapple')
        assert ['apple'] == recipe.missing_structures

    def test_one_missing_structure_long(self):
        recipe = self.parse_command_line('--missing apple')
        assert ['apple'] == recipe.missing_structures

    def test_three_missing_structures(self):
        recipe = self.parse_command_line('--missing apple -mbanana --missing cherry')
        assert ['apple', 'banana', 'cherry'] == recipe.missing_structures

    def test_no_regmat(self):
        recipe = self.parse_minimum_command_line()
        assert not recipe.regmat

    def test_short_regmat(self):
        recipe = self.parse_command_line('-rsomething')
        assert 'something' == recipe.regmat

    def test_long_regmat(self):
        recipe = self.parse_command_line('--regmat something')
        assert 'something' == recipe.regmat

    def test_busy_args(self):
        recipe = self.parse_command_line(
            '--regmat something -iread_me --input read_me_too -owrite_me --threads 7 --missing banana')
        assert ['read_me', 'read_me_too'] == recipe.image_file_names
        assert 'write_me' == recipe.output
        assert 7 == recipe.threads
        assert ['banana'] == recipe.missing_structures
        assert 'something' == recipe.regmat
