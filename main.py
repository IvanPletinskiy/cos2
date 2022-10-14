import math

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

N = 32
KN = lambda N: N // 4
PHASES = {'0': 0, r'\pi/2': math.pi / 16}
EXPECTED_ROOT_MEAN_SQUARE = 0.707
EXPECTED_AMPLITUDE = 1

LABELS = [
    '(1-Фурье)',
    'Погр-СКЗ-1',
    'Погр-СКЗ-2'
]


def calc_fourier_amplitude(sequence):
    N = len(sequence)
    sin_sum, cos_sum = 0, 0
    for i, y in enumerate(sequence):
        angle = 2 * math.pi * i / N
        sin_sum += y * math.sin(angle)
        cos_sum += y * math.cos(angle)
    sin_sum *= 2/N
    cos_sum *= 2/N
    return math.sqrt(sin_sum**2 + cos_sum**2)


def calc_harmonic_value(x, N, phase):
    return math.sin((2 * math.pi * x) / N + phase)


def calc_square_root_mean_1(sequence):
    N = len(sequence)
    return math.sqrt(1 / (N + 1) * sum(x ** 2 for x in sequence))


def calc_square_root_mean_2(sequence):
    N = len(sequence)
    return math.sqrt(1 / (N + 1) * sum(x ** 2 for x in sequence) - (1 / (N + 1) * sum(sequence)) ** 2)


def generate_signal_list(M, N, phase):
    sequence = []
    for i in range(M):
        sequence.append(calc_harmonic_value(i, N, phase))
    return sequence


def calc_parameters(KN, N, phase):
    graphs = {"rmsc_error": [], "rmsd_error": [], "fourier_error": []}
    K = KN(N)
    for x in range(K, 5 * N):
        sequence = generate_signal_list(x, N, phase)

        smr1 = calc_square_root_mean_1(sequence)
        smr2 = calc_square_root_mean_2(sequence)
        fourier = calc_fourier_amplitude(sequence)

        graphs["rmsc_error"].append({"x": x / N, "y": EXPECTED_ROOT_MEAN_SQUARE - smr1})
        graphs["rmsd_error"].append({"x": x / N, "y": EXPECTED_ROOT_MEAN_SQUARE - smr2})
        graphs["fourier_error"].append({"x": x / N, "y": EXPECTED_AMPLITUDE - fourier})

    return graphs


def point_converter(graph):
    xticks = []
    yticks = []
    for point in graph:
        xticks.append(point["x"])
        yticks.append(point["y"])

    return {"x": xticks, "y": yticks}


def graph_controller(K, N):

    def update(x):
        X = int(x)
        l = 0
        for i, phase_unit in enumerate(PHASES.items()):
            desc, phase = phase_unit
            graphs = calc_parameters(KN, X, phase)

            for j, graph in graphs.items():
                points = point_converter(graph)
                plots[l].set_xdata(points["x"])
                plots[l].set_ydata(points["y"])
                l += 1

    plots = []
    fig = plt.figure(figsize=(15, 7))
    plt.grid(True)

    for i, phase_unit in enumerate(PHASES.items()):
        desc, phase = phase_unit
        graphs = calc_parameters(KN, N, phase)
        plt.subplot(1, 2, i + 1)
        l = 0

        for j, graph in graphs.items():
            points = point_converter(graph)
            tempPlot, = plt.plot(points["x"], points["y"], label='${}$'.format(LABELS[l]))
            plots.append(tempPlot)
            l += 1
        plt.legend(loc='best')


    axes = plt.axes([0.37, 0.93, 0.25, 0.03])
    slider = Slider(axes, '${}$'.format('N'), 2, 128, valinit=N, valfmt=r'$%d$')
    slider.on_changed(update)

    return slider

def main():
    keep = graph_controller(KN, N)
    plt.show()


if __name__ == '__main__':
    main()