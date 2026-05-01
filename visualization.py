import numpy as np
import matplotlib.pyplot as plt


def plot_curve(curve, ax=None, x_range=(-3, 4), color="steelblue", label=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    x = np.linspace(x_range[0], x_range[1], 1000)
    rhs = x ** 3 + curve.a * x + curve.b
    valid = rhs >= 0
    y = np.sqrt(np.where(valid, rhs, 0))

    x_valid = x[valid]
    y_valid = y[valid]

    ax.plot(x_valid, y_valid, color=color, label=label)
    ax.plot(x_valid, -y_valid, color=color)

    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return ax


def plot_point_addition(p1, p2, curve, ax=None, x_range=(-3, 4)):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    plot_curve(curve, ax=ax, x_range=x_range)

    x1, y1 = p1
    x2, y2 = p2
    is_doubling = p1 == p2

    if is_doubling:
        slope = (3 * x1 * x1 + curve.a) / (2 * y1)
        title = f"Point Doubling on y² = x³ + {curve.a}x + {curve.b}"
        line_label = "tangent at A"
        result_label = "C = 2A"
    else:
        slope = (y2 - y1) / (x2 - x1)
        title = f"Point Addition on y² = x³ + {curve.a}x + {curve.b}"
        line_label = "line through A and B"
        result_label = "C = A + B"

    x3 = slope * slope - x1 - x2
    neg_y3 = slope * (x3 - x1) + y1
    y3 = -neg_y3

    line_x = np.array([min(x1, x2, x3) - 1, max(x1, x2, x3) + 1])
    line_y = slope * (line_x - x1) + y1
    ax.plot(line_x, line_y, "--", color="orange", label=line_label)

    ax.plot([x3, x3], [neg_y3, y3], ":", color="gray", label="reflection")

    ax.plot(x1, y1, "o", color="red", markersize=10, zorder=5)
    if not is_doubling:
        ax.plot(x2, y2, "o", color="red", markersize=10, zorder=5)
    ax.plot(x3, neg_y3, "o", color="purple", markersize=10, zorder=5)
    ax.plot(x3, y3, "o", color="green", markersize=10, zorder=5)

    ax.annotate("A", (x1, y1), textcoords="offset points", xytext=(8, 8), fontsize=12, fontweight="bold")
    if not is_doubling:
        ax.annotate("B", (x2, y2), textcoords="offset points", xytext=(8, 8), fontsize=12, fontweight="bold")
    ax.annotate("-C", (x3, neg_y3), textcoords="offset points", xytext=(8, -4), fontsize=12, fontweight="bold")
    ax.annotate(result_label, (x3, y3), textcoords="offset points", xytext=(8, 8), fontsize=12, fontweight="bold")

    ax.set_title(title)
    ax.legend(loc="lower right")
    return ax


def plot_curve_points_ff(curve, ax=None, point_color="steelblue", label=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    points = []
    for x in range(curve.p):
        rhs = (x ** 3 + curve.a * x + curve.b) % curve.p
        for y in range(curve.p):
            if (y * y) % curve.p == rhs:
                points.append((x, y))

    xs = [pt[0] for pt in points]
    ys = [pt[1] for pt in points]
    ax.scatter(xs, ys, color=point_color, s=80, zorder=3,
               edgecolors="black", linewidths=0.5, label=label)

    ax.set_xlim(-1, curve.p)
    ax.set_ylim(-1, curve.p)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    step = max(1, curve.p // 17)
    ax.set_xticks(range(0, curve.p, step))
    ax.set_yticks(range(0, curve.p, step))
    return ax


def plot_point_addition_ff(p1, p2, curve, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(8.5, 8.5))

    plot_curve_points_ff(curve, ax=ax, label="curve points")

    p = curve.p
    x1, y1 = p1
    x2, y2 = p2
    is_doubling = p1 == p2

    if is_doubling:
        slope_top = (3 * x1 * x1 + curve.a) % p
        slope_bottom = (2 * y1) % p
        title = f"Point Doubling on y² = x³ + {curve.a}x + {curve.b} (mod {p})"
        line_label = "tangent line"
        result_label = "C = 2A"
    else:
        slope_top = (y2 - y1) % p
        slope_bottom = (x2 - x1) % p
        title = f"Point Addition on y² = x³ + {curve.a}x + {curve.b} (mod {p})"
        line_label = "line through A and B"
        result_label = "C = A + B"

    slope = (slope_top * pow(slope_bottom, -1, p)) % p
    x3 = (slope * slope - x1 - x2) % p
    y3 = (slope * (x1 - x3) - y1) % p
    neg_y3 = (-y3) % p

    # draw the line in real coordinates from A through B to -C, clipped to
    # the field bounds so it stays on the plot
    xs_focus = [x1, x2, x3]
    x_lo = min(xs_focus) - 0.5
    x_hi = max(xs_focus) + 0.5
    line_x = np.linspace(x_lo, x_hi, 300)
    line_y = slope * (line_x - x1) + y1
    in_field = (line_y >= -0.5) & (line_y <= p - 0.5)
    line_y_plot = np.where(in_field, line_y, np.nan)
    ax.plot(line_x, line_y_plot, "--", color="orange", linewidth=2,
            zorder=2, label=line_label)

    # reflection arrow from -C to C
    ax.annotate("", xy=(x3, y3), xytext=(x3, neg_y3),
                arrowprops=dict(arrowstyle="->", color="dimgray",
                                linestyle=":", lw=1.5))

    ax.plot(x1, y1, "o", color="red", markersize=14, zorder=5, markeredgecolor="black")
    if not is_doubling:
        ax.plot(x2, y2, "o", color="red", markersize=14, zorder=5, markeredgecolor="black")
    ax.plot(x3, neg_y3, "o", color="purple", markersize=14, zorder=5, markeredgecolor="black")
    ax.plot(x3, y3, "o", color="green", markersize=14, zorder=5, markeredgecolor="black")

    ax.annotate("A", (x1, y1), textcoords="offset points", xytext=(8, 8),
                fontsize=13, fontweight="bold")
    if not is_doubling:
        ax.annotate("B", (x2, y2), textcoords="offset points", xytext=(8, 8),
                    fontsize=13, fontweight="bold")
    ax.annotate("-C", (x3, neg_y3), textcoords="offset points", xytext=(8, -4),
                fontsize=13, fontweight="bold")
    ax.annotate(result_label, (x3, y3), textcoords="offset points", xytext=(8, 8),
                fontsize=13, fontweight="bold")

    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=9)
    return ax


def plot_scalar_multiplication_ff(generator, k, curve, point_add, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 10))

    plot_curve_points_ff(curve, ax=ax, point_color="lightsteelblue",
                         label="curve points")

    # build the sequence G, 2G, ..., kG by repeatedly adding G
    multiples = []
    P = None
    for _ in range(k):
        P = generator if P is None else point_add(P, generator, curve)
        if P is None:
            break
        multiples.append(P)

    xs = [pt[0] for pt in multiples]
    ys = [pt[1] for pt in multiples]

    ax.plot(xs, ys, "-", color="darkorange", linewidth=0.7, alpha=0.6,
            zorder=2, label=f"path: G → 2G → ... → {len(multiples)}G")

    gx, gy = generator
    kx, ky = multiples[-1]
    ax.plot(gx, gy, "*", color="red", markersize=22, zorder=5,
            markeredgecolor="black", label=f"G = ({gx}, {gy})")
    ax.plot(kx, ky, "o", color="green", markersize=14, zorder=5,
            markeredgecolor="black",
            label=f"{len(multiples)}G = ({kx}, {ky})")

    ax.set_title(
        f"Scalar multiplication on y² = x³ + {curve.a}x + {curve.b}"
        f" (mod {curve.p})"
    )
    return ax
