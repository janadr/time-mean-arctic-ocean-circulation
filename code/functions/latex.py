

def set_size(width, fraction=1, subplots=(1, 1)):
    """
    Set figure dimensions to avoid scaling in LaTeX.
    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 452
    elif width == "AGU":
        width_pt = 397.48499
    elif width == "article":
        width_pt = 221
    elif width == 'beamer':
        width_pt = 307.28987
    elif width == 'poster':
        width_pt = 1192
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

def produce_table(data, hheader=None, vheader=None):
    """
    Spagetthi code producing a vertical laTEX table.
    data has shape of an array
        [[x0, x1, x2, ..., xN],
         [y0, y1, y2, ..., yN],
         ...          ]
    where
    hheader = list/array of horizontal header
    vheader = list/array of vertical header
    If greek letters are used, then they must be enclosed by $$,
    i. e. $\lambda$
    """
    tableString = ""
    n = len(data[:][0])
    tableString += "\\begin{table}[htbp]\n"
    tableString += "\\centering\n"
    tableString += "\\begin{{tabular}}{{l{0:s}}}\n".format("c"*(n-1))
    tableString += "\\hline\n"
    # creating header
    if hheader is not None:
        for element in hheader:
            tableString += f"\\textbf{{{element}}} & "
    tableString = tableString[:-2] + "\\\\\n"
    tableString += "\\hline\n"
    # creating table elements
    for j in range(len(data[0, :])):
        if vheader is not None:
            tableString += f"\\textbf{{{vheader[j]}}} & "
        for i in range(len(data[:, 0])):
            tableString += f"{data[i, j]:.2f} & "
        tableString = tableString[:-2] + "\\\\\n"
    tableString = tableString[:-4] + "\n"
    tableString += "\\end{tabular}\n"
    tableString += "\\caption{}\n"
    tableString += "\\label{table:}\n"
    tableString += "\\end{table}\n"
    return tableString
