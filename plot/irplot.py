#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib as mpl

mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.sans-serif'] = 'Roboto'

def plot_spectrum(ax, filename, color, label, fill=False):
    """Plot a single spectrum from a file."""
    data = np.loadtxt(filename)
    x, y = data[:, 0], data[:, 1]
    print(f"Plotting data from: {filename}")
    
    ax.plot(x, y, label=label, color=color)
    if fill:
        ax.fill_between(x, 0, y, color=color, alpha=0.25)
    return x.min(), x.max(), y.max()

def main():
    parser = argparse.ArgumentParser(description="Plot data from multiple text files.")
    parser.add_argument("files", type=str, nargs='+', help="Data files to plot")
    parser.add_argument("-fill","--fill", action="store_true", help="Fill the area under the plots")
    parser.add_argument("-xmax","--xmax", type=float, help="Set the upper x-axis limit manually")
    
    args = parser.parse_args()
    
    # Colors for each plot
    colors = ['firebrick', 'royalblue', 'green', 'darkorange', 'purple']  # Add more colors as needed
    labels = ['reference', 'computed']  # Extend or modify labels as needed
    
    # Create the figure and axis objects
    fig, ax = plt.figure(figsize=(8, 4.5)), plt.axes()
    fig.subplots_adjust(right=0.97, top=0.99, bottom=0.12)    

    min_x = float('inf')
    max_x = 0
    max_y = 0

    # Plot each file
    for i, filename in enumerate(args.files):
        color = colors[i % len(colors)]  # Cycle through colors if more files than colors
        label = labels[i % len(labels)] if i < len(labels) else f"{filename}"
        x_min, x_max, y_max_current = plot_spectrum(ax, filename, color, label, fill=args.fill)

        # Update x and y limits
        if x_min < min_x:
            min_x = x_min
        if x_max > max_x:
            max_x = x_max
        if y_max_current > max_y:
            max_y = y_max_current
 
    
    # Set axis labels
    ax.set_xlabel('Frequency / cm$^{-1}$')
    ax.set_ylabel('Intensity')
    
    # Set x-axis and y-axis limits
    ax.set_xlim(min_x, args.xmax if args.xmax else max_x)
    ax.set_ylim(0, max_y * 1.05)  # Extend upper y-limit by 5%

    # Add a legend
    ax.legend(fontsize=13)
    
    # Display the plot
    plt.savefig('spectrum.svg')
    plt.savefig('spectrum.pdf',dpi=600)
    print('Saving "spectrum.svg" and "spectrum.pdf"')
    plt.show()

if __name__ == "__main__":
    main()

