from .plot import Plot

class Polygons(Plot):

    def __init__(self, sim_name, name, polygons, linewidth=0.4, color='yellow',
                 **kwargs):
        Plot.__init__(self, sim_name, name, **kwargs)
        self.polygons = polygons
        self.linewidth = linewidth
        self.color = color

    def plot(self, ax, **kwargs):
        for p in self.polygons:
            ax.plot(*p.exterior.xy, linewidth=self.linewidth, color=self.color,
                    **kwargs)


