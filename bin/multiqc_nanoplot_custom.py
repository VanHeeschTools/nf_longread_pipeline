#!/usr/bin/env python

from multiqc.modules.nanoplot import NanoPlotModule
from multiqc.plots import linegraph
import logging
import re
from bs4 import BeautifulSoup
import json

log = logging.getLogger('multiqc')

class CustomNanoPlotModule(NanoPlotModule):
    def __init__(self):
        super(CustomNanoPlotModule, self).__init__()
        self.nanoplot_data = dict()

    def parse_nanoplot_report(self, f):
        super().parse_nanoplot_report(f)
        
        sample = f['s_name']
        content = f['f']
        soup = BeautifulSoup(content, 'html.parser')

        # Extract data for Length vs Quality plot
        quality_data = self.extract_plot_data(soup, 'LengthvsQualityScatterPlot')
        if quality_data:
            self.nanoplot_data[sample] = quality_data

    def extract_plot_data(self, soup, plot_id):
        script_tag = soup.find('script', text=re.compile(plot_id))
        if script_tag:
            match = re.search(r'var data = (\[.*?\]);', script_tag.string, re.DOTALL)
            if match:
                data_str = match.group(1)
                data = json.loads(data_str)
                return {'x': [point['x'] for point in data],
                        'y': [point['y'] for point in data]}
        return None

    def nanoplot_quality_plot(self):
        data = dict()
        for sample, d in self.nanoplot_data.items():
            data[sample] = {
                'x': d['x'],
                'y': d['y']
            }

        return linegraph.plot(data, {
            'id': 'nanoplot_quality_plot',
            'title': 'NanoPlot: Read Length vs Quality',
            'xlab': 'Read Length',
            'ylab': 'Quality score',
            'xLog': True,
            'yMin': 0,
            'tt_label': '<b>{point.x:.0f} bp</b>: {point.y:.2f}',
        })

    def add_quality_plot(self):
        self.add_section(
            name='Read Length vs Quality',
            anchor='nanoplot_quality',
            plot=self.nanoplot_quality_plot()
        )

# This function will be called by MultiQC to initialize the module
def multiqc_nanoplot_custom():
    return CustomNanoPlotModule()