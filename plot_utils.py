import matplotlib.pyplot as plt
import json
import argparse

mom_names = ["intensity", "mean", "sigma", "skewness"]


def plot_moments(prefix):
    data = json.load(open(prefix + ".json"))
    y = data["position"]
    plt.title("Velocity moments")
    for i, m in enumerate(mom_names, start=1):
        plt.subplot(len(mom_names), 1, i)
        plt.plot(y, data[m])
        plt.ylabel(m)
    plt.xlabel("Position along slit, arcsec")
    plt.savefig(prefix + "-plot.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Exercise the plotting routines"""
    )
    parser.add_argument(
        "stampname", type=str,
        help="""Prefix of stamp file (e.g., p84-O_I_5577)"""
    )
    cmd_args = parser.parse_args()
    plot_moments("Stamps/{}-stamp-nc-moments".format(cmd_args.stampname))
