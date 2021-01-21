from __future__ import absolute_import
import matplotlib.pyplot as plt

def set_shared_ylabel(a, ylabel, labelpad = 0.01):
    """Set a y label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    ylabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""

    f = plt.figure()  
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    top = a[0].get_position().y1
    bottom = a[-1].get_position().y0

    # get the coordinates of the left side of the tick labels 
    x0 = 1
    for at in a:
        at.set_ylabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.yaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.inverse_transformed(f.transFigure)
        xt = bboxes.x0
        if xt < x0:
            x0 = xt
    tick_label_left = x0

    # set position of label
    a[-1].set_ylabel(ylabel)
    a[-1].yaxis.set_label_coords(tick_label_left-labelpad, (bottom + top)/2, transform=f.transFigure)

def set_shared_xlabel(a, xlabel, labelpad = 0.01):
    """Set a x label shared by multiple axes
    Parameters
    ----------
    a: list of axes
    ylabel: string
    labelpad: float
        Sets the padding between ticklabels and axis label"""

    f = plt.figure()  
    f.canvas.draw() #sets f.canvas.renderer needed below

    # get the center position for all plots
    left = a[0].get_position().x0
    right = a[-1].get_position().x1

    # get the coordinates of the left side of the tick labels 
    y0 = 1
    for at in a:
        at.set_xlabel('') # just to make sure we don't and up with multiple labels
        bboxes, _ = at.xaxis.get_ticklabel_extents(f.canvas.renderer)
        bboxes = bboxes.inverse_transformed(f.transFigure)
        yt = bboxes.y0
        if yt < y0:
            y0 = yt
    tick_label_bottom = y0

    # set position of label
    a[-1].set_xlabel(xlabel)
    a[-1].xaxis.set_label_coords((left + right)/2, tick_label_bottom-labelpad, transform=f.transFigure)
