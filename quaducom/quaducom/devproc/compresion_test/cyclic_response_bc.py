__author__ = 'campsb'
__date__ = 'Nov. 8, 2017'


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


##########################################################################################################
# variables
##########################################################################################################

test= 'CT_120_1_9'

# directory of script
data_dir= '/users/macretina/sciebo/imb/python'

# directory of raw data
raw_data_dir=os.path.join(data_dir, 'rawdata_CT', '%s.csv' %(test))

# directory of reports
report_dir=os.path.join(data_dir, 'report_CT', '%s' %(test))

# directory of argmax_file
file_dir=os.path.join(report_dir, '%s' %(test))

# creates a directory for report charts and files if there is no such directory
if not os.path.isdir(report_dir):
    os.makedirs(report_dir)


##########################################################################################################
# definition of design of charts
###########################################################################################################

# font and fontsize
font={'name': 'arial','size': 8}

# line widths of main and secondary axes
line_widths = {'main': 1.0, 'sec': 0.5, 'graph':1.0}

# line widths of ticks
ticks={'width_major':1.0, 'width_minor':0.5, 'size_major':3.0, 'size_minor': 3.0, 'log': False}

# size and weight of font of title
font_title={'size': 10, 'weight':'bold'}

# size and weight of font of title
font_axistitle={'size': 8, 'weight':'bold'}

# dimensions (in.) and resolution of chart
fig_dimension={'length':8*0.3937, 'hight':6*0.3937, 'dpi':300}

# position (in.) of chart
adjust_graph={'bottom': 0.21, 'top':0.9, 'right': 0.95, 'left': 0.18}

# font size and spacing of legend
adjust_legend={'size': 6, 'spacing': 0.2}

# domain of x-axis (strain)
adjust_xlim={'xWA_min': 0., 'xWA_max': 1.5}

# domain of y-axis (force)
adjust_ylim={'yF_min': 0.,'yF_max': 1000}


##########################################################################################################
# definition of functions
##########################################################################################################

plt.rc("font", family=font['name'], size=font['size'])

# function to convert force positive
def convert_force_pos(force):
    return force*-1.0

# function to define chart layout
def definiere_design(ax, title, xlabel, ylabel, xlims, ylims, log=ticks['log']):

    # defines chart title
    ax.set_title(title, size=font_title['size'], fontweight=font_title['weight'])

    if log:
        ax.set_xscale('log')

    # defines axis intercepts
    ax.set_xlim(xlims[0], xlims[1])
    ax.set_ylim(ylims[0], ylims[1])

    # label of x-axis
    ax.set_xlabel(xlabel, size=font_axistitle['size'],fontweight=font_axistitle['weight'])

    # label of x-axis
    ax.set_ylabel(ylabel, size=font_axistitle['size'],fontweight=font_axistitle['weight'])

    #orientation of plot
    plt.subplots_adjust(bottom=adjust_graph['bottom'], top=adjust_graph['top'], left=adjust_graph['left'], right=adjust_graph['right'])

    # design of ticks
    ax.tick_params(axis='both', which='major', direction='out', size=ticks['size_major'], width=ticks['width_major'])
    ax.tick_params(axis='both', which='minor', direction='out', size=ticks['size_minor'], width=ticks['width_minor'])

    # position of ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.spines["left"].set_linewidth(line_widths['main'])
    ax.spines["bottom"].set_linewidth(line_widths['main'])
    ax.spines["right"].set_linewidth(line_widths['sec'])
    ax.spines["top"].set_linewidth(line_widths['sec'])
    ax.xaxis.get_major_formatter().set_powerlimits((0,0))


##########################################################################################################
# create data base for chart
##########################################################################################################

# definition of column names
names=['Zeit [s]', 'Kraft F [kN]', 'Maschinenweg [mm]', 'WA_1 [mm]', 'WA_2 [mm]', 'WA_3 [mm]']

# force
F=names[1]
# strain of gauge 1
WA_1=names[3]
# strain of gauge 2
WA_2=names[4]
# strain of gauge 3
WA_3=names[5]

# import of raw data
data=pd.read_csv(raw_data_dir, sep=';', decimal=',', encoding='iso-8859-1', skiprows=2, names=names, nrows=None)


##########################################################################################################
# application of functions
##########################################################################################################

# .map executes defined function fpr every row
data[F]=data[F].map(convert_force_pos)


##########################################################################################################
# supply of plot data
##########################################################################################################

# all data after maximum force will not be considered
data_cut=data[0:data[F].idxmax()+1]

# creates data set with force and strain of gauge 1
WA1_data=data_cut.ix[:,[F, WA_1]]
# deletes empty rows
WA1_data.dropna(axis=0, inplace= True)


WA2_data=data_cut.ix[:,[F, WA_2]]
WA2_data.dropna(axis=0, inplace= True)

WA3_data=data_cut.ix[:,[F, WA_3]]
WA3_data.dropna(axis=0, inplace= True)


##########################################################################################################
# plot and storage of charts
##########################################################################################################

# gauge WA_1

fig=plt.figure(figsize=(fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])

ax1=fig.add_subplot(1,1,1)

WA1_data.plot(x=WA_1, y=F, ax=ax1, color=('0.0'), linewidth=line_widths['graph'], label='WA_1')

definiere_design(ax1, test, 'Verformung w [mm]', 'Kraft F [kN]', (adjust_xlim['xWA_min'], adjust_xlim['xWA_max']),(adjust_ylim['yF_min'], adjust_ylim['yF_max']))

# position legend (1=top right, 2=top left, 3=bottom left, 4=bottom right)
ax1.legend(loc=2, prop={'size':adjust_legend['size']}, labelspacing=adjust_legend['spacing'])

fig.savefig(report_dir + '//WA_1.eps')


# gauge WA_2

fig=plt.figure(figsize=(fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])

ax2=fig.add_subplot(111)

WA2_data.plot(x=WA_2, y=F, ax=ax2, color=('0.0'), linewidth=line_widths['graph'], label='WA_2')

definiere_design(ax2, test, 'Verformung w [mm]', 'Kraft F [kN]', (adjust_xlim['xWA_min'], adjust_xlim['xWA_max']),(adjust_ylim['yF_min'], adjust_ylim['yF_max']))

ax2.legend(loc=2, prop={'size':adjust_legend['size']}, labelspacing=adjust_legend['spacing'])

fig.savefig(report_dir + '//WA_2.eps')


# gauge WA_3

fig=plt.figure(figsize=(fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])

ax3=fig.add_subplot(111)

WA3_data.plot(x=WA_3, y=F, ax=ax3, color=('0.0'), linewidth=line_widths['graph'], label='WA_3')

definiere_design(ax3, test, 'Verformung w [mm]', 'Kraft F [kN]', (adjust_xlim['xWA_min'], adjust_xlim['xWA_max']),(adjust_ylim['yF_min'], adjust_ylim['yF_max']))

ax3.legend(loc=2, prop={'size':adjust_legend['size']}, labelspacing=adjust_legend['spacing'])

fig.savefig(report_dir + '//WA_3.eps')

# gauge WA_1 - WA_3

fig=plt.figure(figsize=(fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])

ax1=fig.add_subplot(111)

WA1_data.plot(x=WA_1, y=F, ax=ax1, color=('0.0'), linewidth=line_widths['graph'], label='WA_1')
WA2_data.plot(x=WA_2, y=F, ax=ax1, color=('0.4'), linewidth=line_widths['graph'], label='WA_2')
WA3_data.plot(x=WA_3, y=F, ax=ax1, color=('0.8'), linewidth=line_widths['graph'], label='WA_3')

definiere_design(ax1, test, 'Verformung w [mm]', 'Kraft F [kN]', (adjust_xlim['xWA_min'], adjust_xlim['xWA_max']),(adjust_ylim['yF_min'], adjust_ylim['yF_max']))

ax1.legend(loc=2, prop={'size':adjust_legend['size']}, labelspacing=adjust_legend['spacing'])

fig.savefig(report_dir + '//WA_alle.eps')

##########################################################################################################
# storage of file with F_max
##########################################################################################################

data_output=data.iloc[data[F].idxmax()]

data_str=str(data_output)

with open(file_dir +'.txt', 'w+') as argmax_file:
    argmax_file.write(data_str)
    argmax_file.close

