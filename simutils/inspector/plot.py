from pathlib import Path


class Plot:

    RESULTS_DIR = Path('results')

    def __init__(self, sim_name, name, extension='pdf'):
        self.sim_name = sim_name
        self.name = name
        self.extension = extension
        self.ensure_directory()

    def rename(self, name):
        self.name = name
        return self

    @property
    def directory(self):
        return self.RESULTS_DIR / self.sim_name / 'plots'

    @property
    def filepath(self):
        return self.directory / f'{self.name}.{self.extension}'

    def ensure_directory(self):
        self.directory.mkdir(parents=True, exist_ok=True)
