""" This script will be for the processing of any data for the 'phase 0'
publication of what I believe is the cbass_classic data set.
"""

import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
# NB the pip cartopy install seems to be broken as it doesn't install the required libararies.
# The solution was to install using conda. conda install cartopy.
# I then had to downgrade shapely to 1.5.17. pip install shapely==1.5.17
from cartopy.mpl.gridliner import Gridliner
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy
import matplotlib.gridspec as gridspec
from matplotlib import collections, patches
import sys
import random
from matplotlib.patches import Rectangle, Polygon, Arrow
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import numpy as np
import re

class SampleOrdinationFigure:
    def __init__(self, dynamic_profile_colour=False):
        self.root_dir = os.path.dirname(os.path.realpath(__file__))
        self.input_base_dir = os.path.join(self.root_dir, 'input')

        # sequences abundances
        self.seq_abund_relative_path = os.path.join(self.input_base_dir, '52_DBV_21022019_2019-06-10_07-15-02.209061.seqs.relative.txt')
        self.seq_rel_df = self._make_seq_rel_df()
        self.ordered_seq_names = self._get_ordered_seq_names()
        self.seq_color_dict = self._set_seq_colour_dict()

        # self.seq_color_dict = self._get_seq_color_dict()

        # profile abundances
        self.profile_abund_relative_path = os.path.join(self.input_base_dir, '52_DBV_21022019_2019-06-10_07-15-02.209061.profiles.relative.txt')
        self.profile_uid_to_profile_name_dict = {}
        self.profile_name_to_profile_uid_dict = {}
        self.prof_rel_df = self._make_prof_rel_df()
        self.ordered_prof_names = self._get_ordered_prof_names()
        if dynamic_profile_colour:
            self.prof_pal = self._get_prof_pal()
            self.prof_color_dict = self._get_prof_color_dict()
        else:
            self.prof_color_dict = self._get_hardcoded_profile_colours()

        # Btwn sample distances
        self.distance_base_dir_samples = os.path.join(self.input_base_dir, 'between_sample_distances')
        self.dist_path_samples = os.path.join(self.distance_base_dir_samples, 'A','2019-06-10_07-15-02.209061.bray_curtis_sample_distances_A.dist')
        self.pcoa_samples_path = os.path.join(self.distance_base_dir_samples, 'A','2019-06-10_07-15-02.209061.bray_curtis_samples_PCoA_coords_A.csv')
        self.sample_uid_to_sample_name_dict = None
        self.sample_name_to_sample_uid_dict = None
        self.sample_dist_df = self._make_sample_dist_df()
        self.sample_pcoa_df = self._make_sample_pcoa_df()

        # Btwn type distances
        self.distance_base_dir_types = os.path.join(self.input_base_dir, 'between_profile_distances', 'A')
        self.dist_path_types = os.path.join(self.distance_base_dir_types,
                                              '2019-06-10_07-15-02.209061.bray_curtis_within_clade_profile_distances_A.dist')
        self.pcoa_types_path = os.path.join(self.distance_base_dir_types, '2019-06-10_07-15-02.209061.bray_curtis_profiles_PCoA_coords_A.csv')
        self.type_uid_to_type_name_dict = None
        self.type_name_to_type_uid_dict = None
        self.type_dist_df = self._make_type_dist_df()
        self.type_pcoa_df = self._make_type_pcoa_df()

        # metainfo
        self.meta_path = os.path.join(self.input_base_dir, 'meta_info.xlsx')
        self.meta_df = self._make_meta_df()

        self.fig_out_path = os.path.join(self.root_dir, 'figures')
        os.makedirs(self.fig_out_path, exist_ok=True)
        # self.gis_input_base_path = os.path.join(self.input_base_dir, 'gis')

        # Figure setup
        # have it set up so that we can transfer this over to the cbass_84 script.
        # we will use GridSpecFromSubplotSpec to set up subplots for each of the sample distances and type distances
        # and then below we will have the seq data. for the main plots we will go 5 down 2 across

        self.fig = plt.figure(figsize=(8, 7))

        gs_rows = 6
        gs_cols = 2
        self.gs = gridspec.GridSpec(gs_rows, gs_cols)

        # sample_distances_ordinations and sample_profile_ordinations
        # four columns for each of the ordinations and one column between each of the ordinations = 19 columns
        # the column in between allows space for the component label.
        self.sample_ordination_sub_gs = gridspec.GridSpecFromSubplotSpec(1, 19, subplot_spec=self.gs[:2, :])
        self.pc1_pc2_sample_dist_ax = plt.subplot(self.sample_ordination_sub_gs[0, :4])
        self.pc1_pc3_sample_dist_ax = plt.subplot(self.sample_ordination_sub_gs[0, 5:9])
        self.pc1_pc2_profile_dist_ax = plt.subplot(self.sample_ordination_sub_gs[0, 10:14])
        self.pc1_pc3_profile_dist_ax = plt.subplot(self.sample_ordination_sub_gs[0, 15:19])

        # sequencing and profile info
        # An area that covers rows 3-5 of the super plot gs, and is split into 5 rows and 19 cols
        self.sample_seq_info_sub_gs = gridspec.GridSpecFromSubplotSpec(6, 19, subplot_spec=self.gs[2:7, :])
        self.sample_seq_info_axarr = [plt.subplot(self.sample_seq_info_sub_gs[:3, :4]),
                                      plt.subplot(self.sample_seq_info_sub_gs[:3, 5:9]),
                                      plt.subplot(self.sample_seq_info_sub_gs[:3, 10:14]),
                                      plt.subplot(self.sample_seq_info_sub_gs[:3, 15:19])]
        self.sample_seq_info_legend_axarr = [
            plt.subplot(self.sample_seq_info_sub_gs[4:5, :9]),
            plt.subplot(self.sample_seq_info_sub_gs[4:5, 10:19])]

        self.sites = ['exposed', 'protected']
        self.experiments = ['CBASS', 'Classic']
        self.sites_location_dict = {'exposed': (39.04470, 22.270300), 'protected':(39.051982,22.26900)}
        self.site_color_dict = {'exposed': 'white', 'protected':'black'}

        self.site_marker_dict = {'protected': 'o', 'exposed': 's'}
        self.classic_cbass_colour_dict = {'Classic': '#808080', 'CBASS': '#FFFFFF'}

    def _set_seq_colour_dict(self):
        """Create the colour dictionary that will be used for plotting by assigning a colour from the colour_palette
        to the most abundant seqs first and after that cycle through the grey_pallette assigning colours
        If we are only going to have a legend that is cols x rows as shown below, then we should only use
        that many colours in the plotting."""
        temp_colour_dict = {}
        predefined_colour_dict = self._get_pre_def_colour_dict()

        for seq_name, color_hash in predefined_colour_dict.items():
            if seq_name in self.ordered_seq_names:
                temp_colour_dict[seq_name] = color_hash

        colour_palette, grey_palette = self._get_colour_lists()

        remaining_seqs = [seq for seq in self.ordered_seq_names if seq not in predefined_colour_dict.keys()]

        for i, seq_name in enumerate(remaining_seqs):
            if i < len(colour_palette):
                temp_colour_dict[seq_name] = colour_palette[i]
            else:
                grey_index = i % len(grey_palette)
                temp_colour_dict[seq_name] = grey_palette[grey_index]

        return temp_colour_dict

    def _get_colour_lists(self):
        colour_palette = self._get_colour_list()
        grey_palette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
        return colour_palette, grey_palette

    @staticmethod
    def _get_pre_def_colour_dict():
        """These are the top 40 most abundnant named sequences. I have hardcoded their color."""
        return {
            'A1': "#FFFF00", 'C3': "#1CE6FF", 'C15': "#FF34FF", 'A1bo': "#FF4A46", 'D1': "#008941",
            'C1': "#006FA6", 'C27': "#A30059", 'D4': "#FFDBE5", 'C3u': "#7A4900", 'C42.2': "#0000A6",
            'A1bp': "#63FFAC", 'C115': "#B79762", 'C1b': "#004D43", 'C1d': "#8FB0FF", 'A1c': "#997D87",
            'C66': "#5A0007", 'A1j': "#809693", 'B1': "#FEFFE6", 'A1k': "#1B4400", 'A4': "#4FC601",
            'A1h': "#3B5DFF", 'C50a': "#4A3B53", 'C39': "#FF2F80", 'C3dc': "#61615A", 'D4c': "#BA0900",
            'C3z': "#6B7900", 'C21': "#00C2A0", 'C116': "#FFAA92", 'A1cc': "#FF90C9", 'C72': "#B903AA",
            'C15cl': "#D16100", 'C31': "#DDEFFF", 'C15cw': "#000035", 'A1bv': "#7B4F4B", 'D6': "#A1C299",
            'A4m': "#300018", 'C42a': "#0AA6D8", 'C15cr': "#013349", 'C50l': "#00846F", 'C42g': "#372101"}

    def _get_hardcoded_profile_colours(self):
        return {'A1g/A1-A1l-A1cr-A1o-A1dp-A1p-A1dq-A1dn': '#ed8ed8', 'A1-A1ds-A1z-A1dr-A1bh': '#8ba4d1', 'A1-A1dm': '#aef686', 'A1-A1z-A1do': '#badffb', 'A1-A1du-A1z-A1ds-A1dr': '#fcad97', 'A1-A1z': '#A48684'}

    def plot_ordination_figure(self):

        color_list, marker_list, x_list = self._plot_pc1_pc2_sample_dists()

        self._plot_pc1_pc3_sample_dists(color_list, marker_list, x_list)

        self._plot_type_dists(ax=self.pc1_pc2_profile_dist_ax, second_pc_label='PC2', second_pc_var_explained='13.9')
        self._plot_type_dists(ax=self.pc1_pc3_profile_dist_ax, second_pc_label='PC3', second_pc_var_explained='00.7')

        self._seq_and_type_plotting_site_ordered()
        plt.tight_layout()
        print('saving .png')
        plt.savefig(os.path.join(self.fig_out_path, 'classic_sample_profile_dists_and_seq_info.png'), dpi=1200)
        print('saving .svg')
        plt.savefig(os.path.join(self.fig_out_path, 'classic_sample_profile_dists_and_seq_info.svg'), dpi=1200)

    def _seq_and_type_plotting_site_ordered(self):
        sample_order = self._get_sample_order()
        self._plot_seq_and_type_ax_site_ordered(sample_order)
        self._plot_seq_and_type_legend(ax=self.sample_seq_info_legend_axarr[0])
        self._plot_seq_and_type_legend(ax=self.sample_seq_info_legend_axarr[1], type_plotting=True)

    def _get_sample_order(self):
        sample_plotting_order_matrix = [[[] for _ in range(len(self.sites))] for _ in
                                        range(len(self.ordered_prof_names))]
        # We will order the samples by mostabund type profile and then by site
        # ugly but fastest is just to go through the df multiple times
        for sample_uid in self.prof_rel_df.index.values.tolist():
            prof_name_ind = self.ordered_prof_names.index(self.prof_rel_df.loc[sample_uid].idxmax())
            site_ind = self.sites.index(self.meta_df.at[self.sample_uid_to_sample_name_dict[sample_uid], 'site'])
            sample_plotting_order_matrix[prof_name_ind][site_ind].append(sample_uid)
        sample_order = []
        for i in range(len(sample_plotting_order_matrix)): # for each profile
            for j in range(len(sample_plotting_order_matrix[i])): # for each site
                sample_order.extend(sample_plotting_order_matrix[i][j])
        return sample_order

    def _plot_seq_and_type_ax_site_ordered(self, sample_order):
        # We plot the first 55 because this is clean break in the ITS2 type profiles
        for i, site in enumerate(self.sites):
            for j, experiment in enumerate(self.experiments):
            # ordered sample list should be the samples of the site in the order of sample_order
                ordered_sample_list = [sample_uid for sample_uid in sample_order if self.meta_df.at[self.sample_uid_to_sample_name_dict[sample_uid], 'site'] == site if self.meta_df.at[self.sample_uid_to_sample_name_dict[sample_uid], 'cbass_classic'] == experiment]
                # Now that we are going to have the samples labelled we need to also sort
                # according to the sample number and then temperature.
                ordered_sample_names = [self.sample_uid_to_sample_name_dict[sample_uid] for sample_uid in ordered_sample_list]
                # names are in the following format
                # S2_E_ST_30B where S2 is the coral sample, E is exposed, ST is short term and 30B is the temperature
                new_sample_names_sorted = []
                for sample_num in range(1,8,1):
                    temp_sub_unsorted = []
                    for sample_name in ordered_sample_names:
                        if f'S{sample_num}' in sample_name:
                            temp_sub_unsorted.append(sample_name)
                    temp_sub_sorted = []
                    for temp in [30,31,33,35,36,39]:
                        for sample_name in temp_sub_unsorted:
                            if str(temp) in sample_name:
                                temp_sub_sorted.append(sample_name)
                    new_sample_names_sorted.extend(temp_sub_sorted)
                # now we need to convert this list of sorted sample names back to a list of UIDs
                sample_name_to_sample_uid_dict = {v:k for k, v in self.sample_uid_to_sample_name_dict.items()}
                ordered_sample_list = [sample_name_to_sample_uid_dict[sample_name] for sample_name in new_sample_names_sorted]
                # print out the sample names in each of the four plots so that carol can annotate the svg
                print(f'{site}_{experiment}')

                num_sampls_first_plot = len(ordered_sample_list)
                width = 1 / num_sampls_first_plot
                if i == 0 and j == 0:
                    y_lab=True
                else:
                    y_lab = False
                self._plot_seqs_on_ax(
                    ordered_sample_list=ordered_sample_list,
                    ax=self.sample_seq_info_axarr[(2*i)+j],
                    width=width,
                    x_ind_list=[_ * width for _ in range(num_sampls_first_plot)],
                    num_samples_in_first_plot=num_sampls_first_plot, y_lab=y_lab, title=f'{site}_{experiment}')

    def _plot_seq_and_type_legend(self, ax, type_plotting=False):
        leg_plotter = LegendPlotter(parent_plotter=self, ax=ax, type_plotting=type_plotting)
        leg_plotter.plot_legend_seqs()

    def format_seq_type_axis(self, ax, num_samples_in_first_plot, ordered_sample_list, width, x_ind_list, y_lab):
        ax.set_xlim(0, len(ordered_sample_list) / num_samples_in_first_plot)
        prof_depth = 0.3
        ax.set_ylim(-(prof_depth + 0.25), 1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([x_left + 0.5*width for x_left in x_ind_list])
        # Chris wants the names in the format of Col1 36 instead of e.g. S1_E_ST_36
        # We will implement this here.
        ax.set_xticklabels([self._convert_sample_name_to_new_format(self.sample_uid_to_sample_name_dict[sample_uid]) for sample_uid in ordered_sample_list], rotation=90, size=4)
        ax.set_yticks([])
        ax.plot([x_ind_list[0], x_ind_list[-1] + width], [0, 0], 'k-', linewidth=ax.spines['right']._linewidth)
        if y_lab:
            ax.text(x=-0.05, y=0.5, s='ITS2 sequences', horizontalalignment='center', verticalalignment='center',
                    rotation='vertical', fontsize='xx-small')
            ax.text(x=-0.08, y=-(prof_depth + 0.2) / 2, s='predicted\nprofile', horizontalalignment='center',
                    verticalalignment='center', rotation='vertical', fontsize='xx-small')
        return prof_depth

    def _convert_sample_name_to_new_format(self, sample_name):
        # We need to extract the number and the temperature from the sample name
        match_object = re.match(pattern=r'S(\d)_\w_\w{2}_{1,2}(\d\d)\w', string=sample_name)
        if match_object:
            return f'Col{match_object.group(1)} {match_object.group(2)}'
        else:
            foo = 'bar'

    def _plot_seqs_on_ax(self, ordered_sample_list, ax, width, x_ind_list, num_samples_in_first_plot, y_lab=True, title=None):
        prof_depth = self.format_seq_type_axis(ax, num_samples_in_first_plot, ordered_sample_list, width, x_ind_list, y_lab)

        ax.set_title(title, fontsize='x-small')
        patches_list = []
        color_list = []


        for sample_uid, x_ind in zip(ordered_sample_list, x_ind_list):
            # for each sample we will start at 0 for the y and then add the height of each bar to this
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            # make the sequence rectangles
            non_zero_sample_series, sample_total = self._get_sample_seq_info(sample_uid)

            self._add_seq_rects_to_patches_list(color_list, non_zero_sample_series, patches_list, sample_total, width,
                                                x_ind)

            # make the profile rectangle
            non_zero_sample_series, sample_total = self._get_sample_type_info(non_zero_sample_series, sample_total,
                                                                              sample_uid)
            self._add_type_rects_to_patches_list(color_list, non_zero_sample_series, patches_list, prof_depth,
                                                 sample_total, width, x_ind)


        self._add_seq_and_type_rects_to_axis(ax, color_list, patches_list)
        # ax.autoscale_view()
        # ax.figure.canvas.draw()

    def _add_seq_and_type_rects_to_axis(self, ax, color_list, patches_list):
        listed_colour_map = ListedColormap(color_list)
        patches_collection = PatchCollection(patches_list, cmap=listed_colour_map)
        patches_collection.set_array(np.arange(len(patches_list)))
        ax.add_collection(patches_collection)

    def _add_type_rects_to_patches_list(self, color_list, non_zero_sample_series, patches_list, prof_depth,
                                        sample_total, width, x_ind):
        bottom = 0
        for ser_index, rel_abund in non_zero_sample_series.iteritems():
            height = (rel_abund / sample_total) * -prof_depth
            patches_list.append(Rectangle(
                (x_ind, bottom), width,
                height, color=self.prof_color_dict[ser_index]))
            bottom += height
            color_list.append(self.prof_color_dict[ser_index])

    def _add_seq_rects_to_patches_list(self, color_list, non_zero_sample_series, patches_list, sample_total, width,
                                       x_ind):
        bottom = 0
        for ser_index, rel_abund in non_zero_sample_series.iteritems():
            height = rel_abund / sample_total
            patches_list.append(Rectangle(
                (x_ind, bottom), width,
                height, color=self.seq_color_dict[ser_index]))
            bottom += height
            color_list.append(self.seq_color_dict[ser_index])

    def _get_sample_type_info(self, non_zero_sample_series, sample_total, sample_uid):
        current_sample_series = self.prof_rel_df.loc[sample_uid]
        non_zero_indices = current_sample_series.nonzero()[0]
        non_zero_sample_series = current_sample_series.iloc[non_zero_indices]
        sample_total = non_zero_sample_series.sum()
        return non_zero_sample_series, sample_total

    def _get_sample_seq_info(self, sample_uid):
        current_sample_series = self.seq_rel_df.loc[sample_uid]
        non_zero_indices = current_sample_series.nonzero()[0]
        non_zero_sample_series = current_sample_series.iloc[non_zero_indices]
        sample_total = non_zero_sample_series.sum()
        return non_zero_sample_series, sample_total

    def _plot_type_dists(self, ax, second_pc_label, second_pc_var_explained):
        x_list = []
        y_list = []
        color_list = []
        marker_list = []
        for i, type_uid in enumerate(self.type_pcoa_df.index.values.tolist()):
            x_list.append(self.type_pcoa_df['PC1'][type_uid])
            y_list.append(self.type_pcoa_df[second_pc_label][type_uid])
            type_name = self.type_uid_to_type_name_dict[type_uid]
            type_color = self.prof_color_dict[type_name]
            color_list.append(type_color)
            marker_list.append('o')
        for x, y, c, m in zip(x_list, y_list, color_list, marker_list):
            ax.scatter(x, y, c=c, marker=m, edgecolors='black', linewidth=0.4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('PC1 - 84.9%')
        ax.set_ylabel(f'{second_pc_label} - {second_pc_var_explained}%')
        return

    def _plot_pc1_pc2_sample_dists(self):
        x_list = []
        y_list = []
        color_list = []
        marker_list = []
        for i, sample_uid in enumerate(self.sample_pcoa_df.index.values.tolist()):
            sample_name = self.sample_uid_to_sample_name_dict[sample_uid]
            classic_cbass = self.meta_df.loc[sample_name, 'cbass_classic']
            if classic_cbass == 'n/a':
                continue
            site = self.meta_df.loc[sample_name, 'site']
            site_color = self.classic_cbass_colour_dict[classic_cbass]
            color_list.append(site_color)
            marker_list.append(self.site_marker_dict[site])
            x_list.append(self.sample_pcoa_df['PC1'][sample_uid])
            y_list.append(self.sample_pcoa_df['PC2'][sample_uid])

        for x, y, c, m in zip(x_list, y_list, color_list, marker_list):
            self.pc1_pc2_sample_dist_ax.scatter(x, y, c=c, marker=m, edgecolors='black', linewidth=0.4, s=20)
        self.pc1_pc2_sample_dist_ax.set_xticks([])
        self.pc1_pc2_sample_dist_ax.set_yticks([])
        self.pc1_pc2_sample_dist_ax.set_xlabel('PC1 - 84.9%')
        self.pc1_pc2_sample_dist_ax.set_ylabel('PC2 - 12.3%')
        return color_list, marker_list, x_list

    def _plot_pc1_pc3_sample_dists(self, color_list, marker_list, x_list):
        y_list = []

        for i, sample_uid in enumerate(self.sample_pcoa_df.index.values.tolist()):
            sample_name = self.sample_uid_to_sample_name_dict[sample_uid]
            classic_cbass = self.meta_df.loc[sample_name, 'cbass_classic']
            if classic_cbass == 'n/a':
                continue
            y_list.append(self.sample_pcoa_df['PC3'][sample_uid])
        for x, y, c, m in zip(x_list, y_list, color_list, marker_list):
            self.pc1_pc3_sample_dist_ax.scatter(x, y, c=c, marker=m, edgecolors='black', linewidth=0.4, s=20)
        self.pc1_pc3_sample_dist_ax.set_xticks([])
        self.pc1_pc3_sample_dist_ax.set_yticks([])
        self.pc1_pc3_sample_dist_ax.set_xlabel('PC1 - 84.9%')
        self.pc1_pc3_sample_dist_ax.set_ylabel('PC3 - 01.1%')

    def _get_prof_pal(self):
        return ['#%02x%02x%02x' % rgb_tup for rgb_tup in
                         self.create_colour_list(mix_col=(255, 255, 255), sq_dist_cutoff=5000, num_cols=6,
                                                 time_out_iterations=100000)]

    def _get_prof_color_dict(self):
        temp_col_dict = {}
        for i, prof_name in enumerate(self.ordered_prof_names):
            temp_col_dict[prof_name] = self.prof_pal[i]
        return temp_col_dict

    def _get_ordered_prof_names(self):
        return self.prof_rel_df.sum().sort_values(ascending=False).index.values.tolist()

    def _get_ordered_seq_names(self):
        return self.seq_rel_df.sum().sort_values(ascending=False).index.values.tolist()

    def _make_prof_rel_df(self):
        with open(self.profile_abund_relative_path, 'r') as f:
            prof_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f]]

        df = pd.DataFrame(prof_data)
        self.profile_uid_to_profile_name_dict = {uid:name for uid, name in zip(df.iloc[0,2:].values.tolist(), df.iloc[6, 2:].values.tolist())}
        self.profile_name_to_profile_uid_dict = {name:uid for uid, name in self.profile_uid_to_profile_name_dict.items()}
        df.drop(index=list(range(6)), inplace=True)
        df.drop(columns=1, inplace=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:-6,:]
        df.set_index('ITS2 type profile', drop=True, inplace=True)
        df.index = df.index.astype('int')
        return df.astype('float')

    def _make_seq_rel_df(self):
        with open(self.seq_abund_relative_path, 'r') as f:
            seq_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f]]

        df = pd.DataFrame(seq_data)
        df.iat[0,0] = 'sample_uid'
        df.columns = df.iloc[0]
        df.drop(index=0, inplace=True)
        df.drop(columns='sample_name', inplace=True)
        df.set_index('sample_uid', drop=True, inplace=True)
        df = df[:-3]
        df.index = df.index.astype('int')
        # Get rid of all of the superflous columns only leaving the seq rel counts
        df = df.iloc[:,20:]

        return df.astype('float')

    def _make_sample_dist_df(self):
        with open(self.dist_path_samples, 'r') as f:
            dist_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f][1:]]

        df = pd.DataFrame(dist_data)
        self.sample_uid_to_sample_name_dict = {int(uid):name for name, uid in zip(df[0].values.tolist(), df[1].values.tolist())}
        self.sample_name_to_sample_uid_dict = {name:uid for uid, name in self.sample_uid_to_sample_name_dict.items()}
        df.drop(columns=0, inplace=True)
        df.set_index(keys=1, drop=True, inplace=True)
        df.index = df.index.astype('int')
        df.columns = df.index.values.tolist()

        return df.astype(dtype='float')

    def _make_sample_pcoa_df(self):
        df = pd.read_csv(filepath_or_buffer=self.pcoa_samples_path, sep=',', header=0, index_col=1, skipfooter=1)
        df.drop(columns='sample', inplace=True)
        df.index = df.index.astype('int')
        return df

    def _make_type_dist_df(self):
        with open(self.dist_path_types, 'r') as f:
            dist_data = [out_line.split('\t') for out_line in [line.rstrip() for line in f][1:]]

        df = pd.DataFrame(dist_data)
        self.type_uid_to_type_name_dict = {int(uid): name for name, uid in
                                               zip(df[0].values.tolist(), df[1].values.tolist())}
        self.type_name_to_type_uid_dict = {name: uid for uid, name in self.type_uid_to_type_name_dict.items()}
        df.drop(columns=0, inplace=True)
        df.set_index(keys=1, drop=True, inplace=True)
        df.index = df.index.astype('int')
        df.columns = df.index.values.tolist()

        return df.astype(dtype='float')

    def _make_type_pcoa_df(self):
        df = pd.read_csv(filepath_or_buffer=self.pcoa_types_path, sep=',', header=0, index_col=1, skipfooter=1)
        df.drop(columns='sample', inplace=True)
        df.index = df.index.astype('int')
        return df

    def _make_meta_df(self):
        df = pd.read_excel(self.meta_path, sheet_name=0, header=1, index_col=False, drop_row=0 )
        df.set_index('sample_name', drop=True, inplace=True)
        # make sure that each of the sample names in the distance file are found in the meta file
        sample_names_in_meta = df.index.values.tolist()
        ordered_sample_tups_list = []
        for sample_name, sample_uid in self.sample_name_to_sample_uid_dict.items():
            if sample_name not in sample_names_in_meta:
                print(f'Sample name: {sample_name} not found')
            else:
                ordered_sample_tups_list.append((sample_name, sample_uid))
        wanted_sample_names = [tup[0] for tup in ordered_sample_tups_list]
        # fix the bad coordinate system and then convert the coordinate to float for theses columns
        for i, sample_name in enumerate(df.index.values.tolist()): # for each row
            current_val_lat = df['collection_latitude'][sample_name]
            if not isinstance(current_val_lat, float):
                if 'N' in current_val_lat:
                    df.iat[i, 9] = float(current_val_lat[2:-1])
                    current_val_lon = df['collection_longitude'][sample_name]
                    df.iat[i, 10] = float(current_val_lon[2:-1])
        df = df.loc[wanted_sample_names,:]
        # make a new columns that is site name
        site_name = []
        for i, sample_name in enumerate(df.index.values.tolist()):
            if df['collection_latitude'][sample_name] == 29.514673:
                site_name.append('eilat')
            elif df['collection_latitude'][sample_name] == 22.26302:
                site_name.append('protected')
            elif df['collection_latitude'][sample_name] == 22.26189:
                site_name.append('exposed')
            elif df['collection_latitude'][sample_name] == 22.303411:
                site_name.append('kaust')
            else:
                sys.exit('Mismatch in latitude')
        df['site'] = site_name
        # make a new column that is classic or cbass
        cbass_classic_list = []
        for i, sample_name in enumerate(df.index.values.tolist()):
            if 'ST' in sample_name:
                cbass_classic_list.append('CBASS')
            elif 'LT' in sample_name:
                cbass_classic_list.append('Classic')
            else:
                cbass_classic_list.append('n/a')
        df['cbass_classic'] = cbass_classic_list

        return df

    @staticmethod
    def _get_colour_list():
        colour_list = [
            "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
            "#0CBD66",
            "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
            "#BEC459",
            "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
            "#FF913F",
            "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7",
            "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
            "#201625",
            "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
            "#CB7E98",
            "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
            "#806C66",
            "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5",
            "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
            "#353339",
            "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
            "#001325",
            "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
            "#8D8546",
            "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F",
            "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
            "#A45B02",
            "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
            "#F4D749",
            "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
            "#C6DC99",
            "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A",
            "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
            "#A38469",
            "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
            "#E4FFFC",
            "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
            "#A76F42",
            "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2",
            "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
            "#868E7E",
            "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
            "#00B57F",
            "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
        return colour_list

    def create_colour_list(
            self, sq_dist_cutoff=None, mix_col=None, num_cols=50,
            time_out_iterations=10000, avoid_black_and_white=True):
        new_colours = []
        min_dist = []
        attempt = 0
        while len(new_colours) < num_cols:
            attempt += 1
            # Check to see if we have run out of iteration attempts to find a colour that fits into the colour space
            if attempt > time_out_iterations:
                sys.exit('Colour generation timed out. We have tried {} iterations of colour generation '
                         'and have not been able to find a colour that fits into your defined colour space.\n'
                         'Please lower the number of colours you are trying to find, '
                         'the minimum distance between them, or both.'.format(attempt))
            if mix_col:
                r = int((random.randint(0, 255) + mix_col[0]) / 2)
                g = int((random.randint(0, 255) + mix_col[1]) / 2)
                b = int((random.randint(0, 255) + mix_col[2]) / 2)
            else:
                r = random.randint(0, 255)
                g = random.randint(0, 255)
                b = random.randint(0, 255)

            # now check to see whether the new colour is within a given distance
            # if the avoids are true also
            good_dist = True
            if sq_dist_cutoff:
                dist_list = []
                for i in range(len(new_colours)):
                    distance = (new_colours[i][0] - r) ** 2 + (new_colours[i][1] - g) ** 2 + (
                                new_colours[i][2] - b) ** 2
                    dist_list.append(distance)
                    if distance < sq_dist_cutoff:
                        good_dist = False
                        break
                # now check against black and white
                d_to_black = (r - 0) ** 2 + (g - 0) ** 2 + (b - 0) ** 2
                d_to_white = (r - 255) ** 2 + (g - 255) ** 2 + (b - 255) ** 2
                if avoid_black_and_white:
                    if d_to_black < sq_dist_cutoff or d_to_white < sq_dist_cutoff:
                        good_dist = False
                if dist_list:
                    min_dist.append(min(dist_list))
            if good_dist:
                new_colours.append((r, g, b))
                attempt = 0

        return new_colours

class LegendPlotter:
    """This class can be used by the SeqStackedBarPlotter and the TypeStackedBarPlotter to handle
    the plotting of the legend subplot.
    """
    def __init__(self, parent_plotter, ax, type_plotting=False):
        # whether we are plotting types
        self.type_plotting=type_plotting
        self.parent_plotter = parent_plotter
        self.ax_to_plot_on = ax
        # legend setup parameters
        if type_plotting:
            self.max_n_rows = 3
            self.max_n_cols = 3
        else:
            self.max_n_rows = 3
            self.max_n_cols = 7
        self.num_leg_cells = self.max_n_rows * self.max_n_cols
        self.y_coord_increments = 100 / self.max_n_rows
        self.leg_box_depth = 2 / 3 * self.y_coord_increments

        self.x_coord_increments = 100 / self.max_n_cols
        if type_plotting:
            self.leg_box_width = self.x_coord_increments /5
        else:
            self.leg_box_width = self.x_coord_increments / 3
        self._set_n_rows_and_last_row_len()
        self.column_count = 0

    def plot_legend_seqs(self):
        self._set_ylim_and_x_lim_and_invert_y_axis()

        self._plot_legend_rows()

        self._remove_frames_from_axis()

    def _plot_legend_rows(self):
        if not self.type_plotting:
            sys.stdout.write(
                f'\nGenerating figure legend for {str(self.num_leg_cells)} most common sequences\n')
        else:
            sys.stdout.write(
                f'\nGenerating figure legend for {str(self.num_leg_cells)} most common ITS2 '
                f'type profiles\n')

        for row_increment in range(min(self.n_rows, self.max_n_rows)):

            if self._this_is_last_row_of_legend(row_increment=row_increment):
                for col_increment in range(self.max_n_cols):
                    self._plot_legend_row(row_increment=row_increment, col_increment=col_increment)

                    self.column_count += 1
            else:
                for col_increment in range(self.last_row_len):
                    self._plot_legend_row(row_increment=row_increment, col_increment=col_increment)

                    self.column_count += 1

    def _set_ylim_and_x_lim_and_invert_y_axis(self):
        # Once we know the number of rows, we can also adjust the y axis limits
        self.ax_to_plot_on.set_xlim(0, 100)
        self.ax_to_plot_on.set_ylim(0, ((self.n_rows - 1) * self.y_coord_increments) + self.leg_box_depth)
        self.ax_to_plot_on.invert_yaxis()

    def _set_n_rows_and_last_row_len(self):
        if not self.type_plotting:  # we are plotting sequences
            col_elements_to_plot = len(self.parent_plotter.ordered_seq_names)
        else:  # we are plotting types
            col_elements_to_plot = len(self.parent_plotter.ordered_prof_names)
        if col_elements_to_plot < self.num_leg_cells:
            if col_elements_to_plot % self.max_n_cols != 0:
                self.n_rows = int(col_elements_to_plot / self.max_n_cols) + 1
                self.last_row_len = col_elements_to_plot % self.max_n_cols
            else:
                self.n_rows = int(col_elements_to_plot / self.max_n_cols)
                self.last_row_len = self.max_n_cols
        else:
            self.n_rows = self.max_n_rows
            self.last_row_len = self.max_n_cols

    def _this_is_last_row_of_legend(self, row_increment):
        return (row_increment + 1) != self.n_rows

    def _remove_frames_from_axis(self):
        self.ax_to_plot_on.set_frame_on(False)
        self.ax_to_plot_on.get_xaxis().set_visible(False)
        self.ax_to_plot_on.get_yaxis().set_visible(False)

    def _plot_legend_row(self, row_increment, col_increment):
        leg_box_x, leg_box_y = self._add_legend_rect(col_increment=col_increment, row_increment=row_increment)
        self._add_legend_text(leg_box_x, leg_box_y)

    def _add_legend_text(self, leg_box_x, leg_box_y):
        text_x = leg_box_x + self.leg_box_width + (0.2 * self.leg_box_width)
        text_y = leg_box_y + (0.5 * self.leg_box_depth)
        if not self.type_plotting:
            self.ax_to_plot_on.text(
                text_x, text_y, self.parent_plotter.ordered_seq_names[self.column_count],
                verticalalignment='center', fontsize='xx-small')
        else:
            type_name_to_print = self.parent_plotter.ordered_prof_names[self.column_count]
            if len(type_name_to_print) > 15:
                type_name_to_print = f'{type_name_to_print[:14]}...'

            self.ax_to_plot_on.text(
                text_x, text_y, type_name_to_print,
                verticalalignment='center', fontsize='xx-small')

    def _add_legend_rect(self, col_increment, row_increment):
        leg_box_x = col_increment * self.x_coord_increments
        leg_box_y = row_increment * self.y_coord_increments
        if not self.type_plotting:
            self.ax_to_plot_on.add_patch(Rectangle(
                (leg_box_x, leg_box_y), width=self.leg_box_width, height=self.leg_box_depth,
                color=self.parent_plotter.seq_color_dict[
                    self.parent_plotter.ordered_seq_names[self.column_count]]))
        else:
            self.ax_to_plot_on.add_patch(Rectangle(
                (leg_box_x, leg_box_y), width=self.leg_box_width, height=self.leg_box_depth,
                color=self.parent_plotter.prof_color_dict[
                    self.parent_plotter.ordered_prof_names[self.column_count]]))
        return leg_box_x, leg_box_y

sof = SampleOrdinationFigure()
sof.plot_ordination_figure()
