import ap
import matplotlib.pyplot as plt

def get_times(dict, index):
    result = []
    result_keys = []
    for key, value in dict.items():
        if value.strip():
            print(value)
            result_keys.append(key)
            exercise_average = 0
            exercise_count = 0
            elements = value.split(";")
            for e in elements[1:]:
                print(e)
                exercise_average += float(e.split(",")[index-1])
                exercise_count += 1
            exercise_average /= exercise_count
            result.append(exercise_average)
    return result_keys, result


timing_dictionary = {}
with open("timings_new.txt", "r") as f:
    for line in f:
        (key, value) = line.split(" : ")
        timing_dictionary[key] = value


# Get the data
exercises, graph_info = get_times(timing_dictionary, 1)
exercises, overall_time = get_times(timing_dictionary, 2)
exercises, dg_time = get_times(timing_dictionary, 3)
exercises, fp_time = get_times(timing_dictionary, 4)
exercises, fs_time = get_times(timing_dictionary, 5)


# Create exercise Labels:
labels = [x[-1] for x in exercises]
for i in range(len(labels)):
    label = labels[i] + " " + graph_info[i]

# Create the plots
fig, axs = plt.subplots(2,2)
fig.set_size_inches(18.5, 10.5)

axs[0,0].bar(exercises, overall_time, width=0.75, color="blue")
axs[0,0].set_title('Overall Runtime', fontstyle='italic')
# axs[0,0].set_yscale('log')

axs[0,1].bar(exercises, dg_time, width=0.75, color="purple")
axs[0,1].set_title('Directed Graph Runtime', fontstyle='italic')
# axs[0,1].set_yscale('log')

axs[1,0].bar(exercises, fp_time, width=0.75, color="red")
axs[1,0].set_title('Flow Pathway Runtime', fontstyle='italic')
# axs[1,0].set_yscale('log')

axs[1,1].bar(exercises, fs_time, width=0.75, color="green")
axs[1,1].set_title('Flow Sloution Runtime', fontstyle='italic')
# axs[1,1].set_yscale('log')

plt.show()
# plt.legend()
plt.savefig("results_non_log.pdf")