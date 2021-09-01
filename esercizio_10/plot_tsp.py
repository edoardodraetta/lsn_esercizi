import matplotlib


def Draw_Path(p1, p2, ax):
    x_values = [p1[0], p2[0]]
    y_values = [p1[1], p2[1]]
    ax.plot(x_values, y_values, color="#343f56")


def Plot_TSP(path, X, Y, R, mode, ax=None, label_cities=True):

    if ax is None:
        ax = matplotlib.pyplot.gca()

    if mode == "circle":  # draw the circle
        circle = matplotlib.pyplot.Circle((0, 0), R)
        circle.set_fill(False)
        ax.add_patch(circle)
    if mode == "square":  # draw the square
        square = matplotlib.pyplot.Rectangle((-R, -R), 2 * R, 2 * R, lw=2, color="#2940d3")
        square.set_fill(False)
        ax.add_patch(square)

    # First leg and last leg
    n_cities = len(path)
    p1 = [X[0], Y[0]]
    p2 = [X[path[0]], Y[path[0]]]
    Draw_Path(p1, p2, ax)
    p1 = [X[path[n_cities - 1]], Y[path[n_cities - 1]]]
    p2 = [X[0], Y[0]]
    Draw_Path(p1, p2, ax)

    # rest of the path
    n_cities = len(path)
    for i in range(n_cities - 1):
        a = path[i]
        b = path[i + 1]
        p1 = [X[a], Y[a]]
        p2 = [X[b], Y[b]]
        Draw_Path(p1, p2, ax)

    # plot the cities
    ax.scatter(X[1:], Y[1:])
    ax.scatter(X[0], Y[0], color="#fb9300", alpha=1, label="Starting City")

    if label_cities is True:
        # text placement
        for i in range(len(X) - 1):
            if X[i] > 0:
                dx = R * 0.1
            if X[i] < 0:
                dx = -R * 0.1
            if X[i] > 0:
                dy = R * 0.2
            if X[i] < 0:
                dy = -R * 0.2
            ax.text(X[i] + dx, Y[i] + dy, str(i), color="r")

    # Grid and limits
    ax.grid(True)
    lim = R * 1.5
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])

    # legend
    ax.legend(loc="lower left")

    return ax
