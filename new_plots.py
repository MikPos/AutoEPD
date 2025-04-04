import ap
import matplotlib.pyplot as plt

def read_file_to_dict(file_path):
    data_dict = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            # Assuming each line is formatted as "key: value1,value2,value3,..."
            key, values = line.strip().split(":")
            values_list = values.split(",")  # Split the values by commas into a list
            if len(values_list) == 7:
                data_dict[key] = values_list
            else:
                print(f"Skipping line with incorrect number of values: {line}")
    
    return data_dict

def split_values_to_lists(data_dict):
    lists = [[] for _ in range(7)]  # Create 7 empty lists
    
    for key, values in data_dict.items():
        for i in range(7):
            if values:
                lists[i].append(values[i])  # Append each value to the respective list
    
    return lists

# Create the plots
def plot_lists(lists):
    # Convert string values to floats for plotting (assuming all values are numeric)
    lists = [[float(val) for val in list_] for list_ in lists]

    fig, ax1 = plt.subplots()
    fig.set_figwidth(7)
    fig.set_figheight(4)
    ax1.set_xlabel("Exercise")
    ax1.set_ylabel("Count", color="blue")
    
    # Plot each list with a different color
    ax1.plot(data_dict.keys(), lists[0], label=f"Graph Vertices", color = "blue")
    ax1.plot(data_dict.keys(), lists[1], label=f"Graph Edges", color = "blue", linestyle = 'dotted')
    # ax1.plot(data_dict.keys(), lists[2], label=f"Digraph Vertices", color = "green", linestyle = 'dashed')
    # ax1.plot(data_dict.keys(), lists[3], label=f"Digraph Edges", color = "green", linestyle = 'dashdot')

    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
    ax2.set_ylabel("Seconds", color="red")

    ax2.plot(data_dict.keys(), lists[4], label=f"Overall Runtime", color = "red")
    ax2.plot(data_dict.keys(), lists[5], label=f"Digraph Runtime", color = "red", linestyle = 'dotted')
    ax2.plot(data_dict.keys(), lists[6], label=f"Flow Runtime", color = "red", linestyle = 'dashed')
    
    plt.title('Comparison of Runtime data and graph size.')  # Plot title
    
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=(1.2,0.7))

    plt.grid(True)  # Enable grid for better readability
    plt.tight_layout()
    plt.show()  # Display the plot
    plt.savefig("complexity_data.pdf")


# Example usage:
file_path = 'timings_new.txt'  # Replace with the path to your file
data_dict = read_file_to_dict(file_path)
lists = split_values_to_lists(data_dict)
plot_lists(lists)

# Printing the lists for verification
for i, list_ in enumerate(lists, start=1):
    print(f"List {i}: {list_}")
