'''
Created on 25.10.2013

@author: alexander
'''
def format_plot(axes, xlim=None, ylim=None, xlabel='', ylabel=''):
    '''format 2d-plot black and with with times font 
    '''
    #-------------------------------------------------------------------
    # configure the style of the font to be used for labels and ticks
    #-------------------------------------------------------------------
    #
    from matplotlib.font_manager import FontProperties
    font = FontProperties()
    font.set_name('Script MT')
    font.set_family('serif')
    font.set_style('normal')
#    font.set_size('small')
    font.set_size('large')
    font.set_variant('normal')
    font.set_weight('medium')
    
    if xlim != None and ylim != None:
        axes.axis([0, xlim, 0., ylim], fontproperties=font)

    # format ticks for plot
    #
    locs, labels = axes.xticks()
    axes.xticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.xlabel(xlabel, fontproperties=font)

    locs, labels = axes.yticks()
    axes.yticks(locs, map(lambda x: "%.0f" % x, locs), fontproperties=font)
    axes.ylabel(ylabel, fontproperties=font)