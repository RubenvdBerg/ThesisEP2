import warnings
import matplotlib.pyplot as plt


def only_one_none(a, b, c):
    a, b, c = a is None, b is None, c is None
    return a ^ b ^ c ^ all((a, b, c))


def multi_legend(axes: tuple[plt.Axes,...], **kwargs):
    lines_list = []
    labels_list = []
    for ax in axes:
        lines, labels = ax.get_legend_handles_labels()
        for line, label in zip(lines, labels):
            lines_list.append(line)
            labels_list.append(label)
    axes[0].legend(lines_list, labels_list, **kwargs)
