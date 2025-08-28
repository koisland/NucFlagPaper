import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Generator
from matplotlib.transforms import Transform
from collections import deque


# https://python-graph-gallery.com/591-arrows-with-inflexion-point/
def arrow_inflexion(
    ax: Axes,
    start: tuple[int, int],
    end: tuple[int, int],
    angleA: int,
    angleB: int,
    radius: int = 0,
    color: str = "black",
    transform: Transform | None = None,
    arrow_style: str = "->",
):
    # get the coordinates
    x1, y1 = end
    x2, y2 = start

    # avoid division by zero
    epsilon = 1e-6
    if x2 == x1:
        x2 += epsilon
    if y2 == y1:
        y2 += epsilon

    # select right coordinates
    if transform is None:
        transform = ax.transData

    # add the arrow
    connectionstyle = f"angle,angleA={angleA},angleB={angleB},rad={radius}"
    ax.annotate(
        "",
        xy=(x1, y1),
        xycoords=transform,
        xytext=(x2, y2),
        textcoords=transform,
        size=20,
        arrowprops=dict(
            color=color,
            arrowstyle=arrow_style,
            shrinkA=5,
            shrinkB=5,
            patchA=None,
            patchB=None,
            connectionstyle=connectionstyle,
        ),
    )


def main():
    N_ROW_COL = 14
    fig, axes = plt.subplots(
        nrows=N_ROW_COL,
        ncols=N_ROW_COL,
        figsize=(N_ROW_COL, N_ROW_COL),
        squeeze=True,
        layout="constrained",
    )

    # Identity diagonal
    coord_diag = set(zip(range(N_ROW_COL), range(N_ROW_COL)))

    def tri_coords(
        dim: int, *, offset: int = 0, top: bool = True
    ) -> Generator[tuple[int, int], None, None]:
        for i, x in enumerate(range(dim, 0, -1)):
            for y in range(dim - i):
                if top:
                    yield y + offset, x + offset
                else:
                    yield x + offset, y + offset

    coords_grp_1 = set(tri_coords(dim=3, offset=0))
    coords_grp_2 = set(tri_coords(dim=6, offset=7))
    coords_grps_similar = set((x, y) for y in range(7, N_ROW_COL) for x in range(0, 4))
    coords_traveled = coords_grp_1.union(coords_grp_2)
    data_coords = {}
    for x, y in (*coord_diag, *coords_grp_1, *coords_grp_2, *coords_grps_similar):
        if not data_coords.get(x):
            data_coords[x] = {}
        data_coords[x][y] = 1.0

    # TODO: Draw whole figure first
    # TODO: Then try to make animation (https://matplotlib.org/stable/users/explain/animations/animations.html)
    for x in range(N_ROW_COL):
        for y in range(N_ROW_COL):
            ax: Axes = axes[x, y]
            # Make each ax 3x3
            ax.set_xlim(0, 2)
            ax.set_ylim(0, 2)

            xlabel, ylabel = None, None

            if x == 0 and y == 0:
                xlabel, ylabel = "0", "0"
            elif x == 0:
                xlabel = f"{y}"
            elif y == 0:
                ylabel = f"{x}"

            # Minimize all ticks.
            ax.set_xticks([], [])
            ax.set_yticks([], [])

            if xlabel:
                # Top not bottom.
                ax.xaxis.tick_top()
                ax.set_xlabel(xlabel)
                ax.xaxis.set_label_position("top")
            if ylabel:
                # Horizontal ylabel with padding starting at right of label.
                ax.set_ylabel(ylabel, rotation=0, ha="right")

            if (x, y) in coord_diag:
                ax.set_facecolor("red")
            elif (x, y) in coords_grp_1 or (y, x) in coords_grp_1:
                ax.set_facecolor("#f4ccccff")
            elif (x, y) in coords_grp_2 or (y, x) in coords_grp_2:
                ax.set_facecolor("#c9daf8ff")
            elif (x, y) in coords_grps_similar or (y, x) in coords_grps_similar:
                ax.set_facecolor("#D3D3D3")

    traveled = set()
    window = 5000
    for x in data_coords.keys():
        y = x
        if (x, y) in traveled:
            continue

        positions = deque([(x, y)])
        idents = []
        max_x = x
        while positions:
            nx, ny = positions.popleft()
            if (nx, ny) in traveled:
                continue

            traveled.add((nx, ny))

            if data_coords.get(nx):
                col = data_coords[nx]
            else:
                max_x = nx
                continue

            if col.get(ny):
                ident = col[ny]
            else:
                continue

            try:
                ax: Axes = axes[nx, ny] if nx != ny else None
            except IndexError:
                ax = None
            coord_next_down = (nx, ny + 1)
            coord_next_right = (nx + 1, ny)
            coord_brk_up_left = (nx - 1, ny + 1)
            coord_brk_down_right = (nx + 1, ny - 1)

            if data_coords.get(coord_brk_up_left[0]):
                brk_up_left = data_coords[coord_brk_up_left[0]].get(
                    coord_brk_up_left[1]
                )
            else:
                brk_up_left = None

            if data_coords.get(coord_brk_down_right[0]):
                brk_down_right = data_coords[coord_brk_down_right[0]].get(
                    coord_brk_down_right[1]
                )
            else:
                brk_down_right = None

            positions.appendleft(coord_next_down)
            positions.appendleft(coord_next_right)
            idents.append(ident)

            # Draw arrow.
            # Coords reversed because want on top.
            if (
                ax
                and coord_next_down in coords_traveled
                and coord_next_right in coords_traveled
            ):
                arrow_inflexion(
                    ax=ax,
                    start=(2, 1),
                    end=(1, 0),
                    angleA=0,
                    angleB=90,
                    radius=0,
                    arrow_style="<|-|>",
                )
            elif ax and coord_next_down in coords_traveled:
                arrow_inflexion(
                    ax=ax,
                    start=(1, 1),
                    end=(2, 1),
                    angleA=0,
                    angleB=90,
                    radius=0,
                    arrow_style="-|>",
                )
            elif ax and coord_next_right in coords_traveled:
                arrow_inflexion(
                    ax=ax,
                    start=(1, 1),
                    end=(1, 0),
                    angleA=0,
                    angleB=90,
                    radius=0,
                    arrow_style="-|>",
                )

            if (
                ax
                and not brk_up_left
                and not brk_down_right
                and coord_next_down not in coords_traveled
                and coord_next_right not in coords_traveled
            ):
                max_x = nx
                # Hit end.
                circle = plt.Circle((1, 1), radius=0.25, color="black")
                ax.add_artist(circle)

        if not idents:
            continue

        start = x * window + 1
        end = max_x * window + 1
        length = end - start
        if length <= window:
            continue
        print(start, end)

    fig.savefig("test.png")


if __name__ == "__main__":
    raise SystemExit(main())
