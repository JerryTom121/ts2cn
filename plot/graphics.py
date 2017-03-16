# code to plot graphics

# import libraries
from igraph import Histogram, plot
from numpy import array, median, std, arange
from ts2cn.cn import measures as ms
from matplotlib import pyplot as plt
import imp


def plot_scores(scores, x=None, labels=None, leg_loc='best', x_label='', y_label='', title='', style=None, display='plot',
				filename='', directory='', figsize=()):
	"""
####################################################################################
 Plot several scores together
####################################################################################
scores		List of tuples. Each tuple has several scores to plot
x			List of numbers. The values of the x axis
labels		List of string. The label of each score 
leg_loc		Number or String. Legend location, values 'best' or 0, 'upper right' or 1, 
			'upper left' or 2, 'lower left' or 3, 'lower right' or 4, 'right' or 5, 
			'center left' or 6, 'center right' or 7, 'lower center' or 8, 
			'upper center' or 9, 'center' or 10
x_label		String. X axis label
y_label 	String. Y axis label
tite		String. The plot's title
style		List of tuples. Each tuple has a color, marker and linestyle for each score
			Default values ['r', 'o', '-'], see http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
			for more details
display		String. 'plot' to display the generated plots, 'subplot' for a subplot graphic, something else to save
			to a file
filename	String. The filename with extension
directory	String. The directory where to save the image, remember to end with a slash '/'
figsize		Tuple. The figure size (widht, height) in inches
####################################################################################

usage examples: 
lb='ab'; data=tss[0][4]['F'][:100]; imp.reload(mt); pp=datetime.now(); levels=3; w_size=int(40/2**levels); labels = [[lb]*int(len(data[0])/w_size)] * len(data); wavelet='db4'; print(pp); oo = mt.automated_categorization(np.hstack(baseline), data, labels, w_size, wavelet, levels, seiz_dur, spike_amp, clean_event_window, TFs[:20], silent=True);

dic='/tmp/x/'; figsize=(); titles=['Normal', 'Seizure', 'Spikes', 'Abnormal']; imp.reload(graphics); 
for index in [0,1,2,3]:
	scores=[tuple(j[index] if type(j) == type(np.ndarray([1])) else j for j in i.values()) for i in oo['metrics']]; 
	graphics.plot_scores(scores, TFs, labels=['Accuracy', 'Precision', 'Recall', 'F1'], x_label='Threshold Values', y_label='Score', style=[['r', 'o', '-'], ['b', '*', '-'], ['g', 'v', '-'], ['y', '^', '-']], leg_loc='best', filename=titles[index]+'.png', directory=dic, figsize=figsize, display='save');
	"""
	imp.reload(plt)
	# set the default values
	if not x:
		x = range(len(scores))
	if not labels:
		labels = [''] * len(scores[0])
	if not style:
		style = [['r', 'o', '-']] * len(scores[0])
	
	if figsize:
		plt.figure(figsize=figsize)
	
	for i in range(len(scores[0])):
		plt.plot(x, [j[i] for j in scores], label=labels[i], color=style[i][0], marker=style[i][1], linestyle=style[i][2])
	
	# get an epsilon 1% of the max range
	x_max = max(x)
	x_min = min(x)
	epsilon_x = (x_max - x_min )/ 20
	y_max = max([j for i in scores for j in i])
	y_min = min([j for i in scores for j in i])
	epsilon_y = (y_max - y_min) / 20
	
	plt.xlim(x_min - epsilon_x, x_max + epsilon_x)
	plt.ylim(y_min - epsilon_y, y_max + epsilon_y)
	plt.legend(loc=leg_loc)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.title(title)
	
	if display == 'plot':
		plt.show()
	elif display != 'subplot':
		plt.savefig(directory+filename, format='png')
		plt.close()
	





# TODO include legend of the number of seconds and the amplitude
# TODO allow to set other style parameters other than the color
def plot_sig_pred_grd(signal, sig_xlabel='', sig_ylabel='', sig_style='', grd_truth=[], grd_ylabel='',
				predicted=[], pred_ylabel='', labels={}, labels_style=[], labels_pos=-600, window=1, 
				bar_height=10, v_offset=1, ylim=None, grid_on='on', grid_which='both', grid_axes='x',
				ticks_space=100000, title='', display='plot', filename='', directory='', figsize=()):
	"""
####################################################################################
Plot signal along side the predicted values and groundtruth sample by sample
####################################################################################
signal		List. The signal to plot
sig_xlabel	String. A label for the signal (xlabel)
sig_ylabel	String. A label for the signal (ylabel)
sig_style	String. The style for the signal
grd_truth	List. The list of the ground truth values.
grd_ylabel	String. The label for the ground truth subplot
predicted	List. The list of the predicted values
pred_ylabel	String. The label for the predicted subplot
labels		Dictionary. A dictionary for each predicted or ground truth value 
			to a given label. Both ground truth and predicted values should have
			the same set of labels
labels_pos	Number. The x position of the labels, ground truth and predicted
window		Number. The window size in samples
bar_height	Number. The height of the horizontal bars
v_offset	Number. The vertical offset of the bars
ylim		Tuple. The ylim min and max
grid_on		String. Parameter for function pyplot.grid which sets the grid on or off. Values 'on' (default) or 'off'
grid_which	String. Parameter for function pyplot.grid which sets the type of grid. Values 'major', 'minor', or 'both' (default)
grid_axes	String. Parameter for function pyplot.grid which sets the grid axes. Values 'both', 'x' (default) or 'y'
ticks_space	Number. The space between x ticks 
tite		String. The plot's title
display		String. 'plot' to display the generated plots, 'subplot' for a subplot graphic, something else to save
			to a file
filename	String. The filename with extension
directory	String. The directory where to save the image, remember to end with a slash '/'
figsize		Tuple. The figure size (widht, height) in inches
####################################################################################

Usage example:
graphics.plot_(data[0], 'sig_xlabel', 'sig_label', 'r-', ['ab']*len(oo['events'][0][0]), {'ab': 'Abnormal'}, 'blue', 'Ground Truth',  oo['events'][0][0], {'ab': 'Abnormal', 'no': 'Normal'}, 'green', 'Predicted', {'ab': 'Abnormal', 'no': 'Normal'}, window=5)


imp.reload(graphics); graphics.plot_(data[0], 'sig_xlabel', 'sig_label', 'r-', ['ab']*len(oo['events'][0][0]), {'ab': 'Abnormal'}, ['blue', 'green'], 'Ground Truth',  oo['events'][0][0], {'ab': 'Abnormal', 'no': 'Normal'}, ['blue', 'green'], 'Predicted', {'ab': 'Abnormal', 'no': 'Normal'}, window=5)


lb = 'ab'; signal = data[0]; predicted = oo['events'][0][0]; grd_truth = [lb]*len(predicted); labels_dic = {'ab': 'Abnormal', 'no': 'Normal', 'sp':'Spike', 'sz':'Seizure'}; labels_style = ['red', 'blue', 'green', 'yellow']; imp.reload(graphics); graphics.plot_(signal, 'Samples', 'Signal', 'r-', grd_truth, labels_dic, labels_style, 'Ground Truth',  predicted, labels_dic, labels_style, 'Predicted', labels_dic, window=w_size)

lb = 'ab'; signal = data[0]; predicted = oo['events'][0][0]; grd_truth = [lb]*len(predicted); labels_dic = {'ab': 'Abnormal', 'no': 'Normal', 'sp':'Spike', 'sz':'Seizure'}; labels_style = ['red', 'blue', 'green', 'yellow']; imp.reload(graphics); graphics.plot_(signal, 'Samples', 'Signal', 'r-', grd_truth, labels_dic, labels_style, 'Ground Truth',  predicted, labels_dic, labels_style, 'Predicted', labels_dic, window=w_size, bar_height=5, v_offset=0, display='subplot')


lb = 'ab'; signal = data[0]; predicted = oo['events'][0][0]; grd_truth = [lb]*len(predicted); labels_dic = {'ab': 'Abnormal', 'no': 'Normal', 'sp':'Spike', 'sz':'Seizure'}; labels_style = ['red', 'blue', 'green', 'yellow']; imp.reload(graphics); graphics.plot_(signal, 'Samples', 'Signal', 'c-', grd_truth, 'Ground Truth',  predicted, 'Predicted', labels_dic, labels_pos=-600, window=w_size, bar_height=10, v_offset=10, display='plot')


lb = 'ab'; signal = data[0]; predicted = oo['events'][0][0]; grd_truth = [lb]*len(predicted); labels_dic = {'ab': 'Abnormal', 'no': 'Normal', 'sp':'Spike', 'sz':'Seizure'}; labels_style = ['red', 'blue', 'green', 'yellow']; imp.reload(graphics); graphics.plot_(signal, 'Samples', 'Signal', 'c-', grd_truth, 'Ground Truth',  predicted, 'Predicted', labels_dic, labels_style, labels_pos=-600, window=w_size, bar_height=10, v_offset=0, display='plot', figsize=(15,5))
	"""
	imp.reload(plt)

	import numpy as np
	
	grd_tryth = np.array(grd_truth)
	predicted = np.array(predicted)
	
	if figsize:
		plt.figure(figsize=figsize)
	#plt.subplot(121)
	plt.subplot2grid((1,10), (0,0), colspan=9)
	
	sig_max = 0
	sig_min = 0
	if display != 'subplot':
		sig_max = max(signal)
		sig_min = min(signal)
	
	if ylim:
		sig_min = ylim[0]
		sig_max = ylim[1]
	

	# plot the ground truht horizontal bars
	# ground_truth label
	plt.text(x=labels_pos, y= - v_offset + sig_max, s=20, text=grd_ylabel, rotation='vertical')
	for i in enumerate(labels):
		
		# get the contiguous elements in the grd_truth
		index = np.argwhere(grd_truth == i[1])
		if not index.any(): 
			continue
		size = 1
		pos = index[0]
		positions = []
		for j in range(1, len(index)):
			if index[j-1] +1 == index[j]:
				size += 1
			else:
				plt.barh(-bar_height*i[0] - v_offset + sig_max, width=size*window, left=pos*window, height=bar_height, 
						color=labels_style[i[0]], alpha=0.5, linewidth=0)
				#positions.append((pos, size))
				size = 1
				pos = index[j]
		# add the last one, or the first if there is only one
		#positions.append((pos, size))
		plt.barh(-bar_height*i[0] - v_offset + sig_max, width=size*window, left=pos*window, height=bar_height, 
				color=labels_style[i[0]], alpha=0.5, linewidth=0)
		
		"""
		for sample in [index*window for index,j in enumerate(grd_truth) if j == i[1]]:
			print(sample)
			plt.barh(-bar_height*i[0] - v_offset + sig_max, width=window, left=sample, height=bar_height, 
						color=labels_style[i[0]], alpha=0.5, linewidth=0)
		"""
	# plot the signal
	plt.plot(signal, sig_style, color='cyan', zorder=-1)
	plt.xticks(np.arange(0, len(signal), ticks_space))
	
	from matplotlib.ticker import MultipleLocator
	minorLocator = MultipleLocator(ticks_space/4)
	axes = plt.gca()
	axes.xaxis.set_minor_locator(minorLocator)
	
	plt.grid(b=grid_on, which=grid_which, axis=grid_axes)
#	plt.tick_params(which=grid_which, axis=grid_axes)
	plt.plot([0,len(signal)], [0, 0], color='gray', zorder=-2)
	plt.xlabel(sig_xlabel)
	plt.ylabel(sig_ylabel)
	plt.xlim(0, len(signal))
	if ylim:
		plt.ylim(ylim[0], ylim[1])

	# plot the predicted values
	# predicted values label
	plt.text(x=labels_pos, y= v_offset + sig_min, s='normal', text=pred_ylabel, rotation='vertical')
	

	for i in enumerate(labels):
		

		index = np.argwhere(predicted == i[1])
		if not index.any(): 
			continue
		size = 1
		pos = index[0]
		positions = []
		for j in range(1, len(index)):
			if index[j-1] +1 == index[j]:
				size += 1
			else:
				plt.barh(bar_height*i[0] + v_offset + sig_min, width=size*window, left=pos*window, height=bar_height, 
						color=labels_style[i[0]], alpha=0.5, linewidth=0)
				#positions.append((pos, size))
				size = 1
				pos = index[j]

		# add the last one, or the first if there is only one
		plt.barh(bar_height*i[0] + v_offset + sig_min, width=size*window, left=pos*window, height=bar_height, 
				color=labels_style[i[0]], alpha=0.5, linewidth=0)
		#positions.append((pos, size))
		"""
		for sample in [index*window for index,j in enumerate(predicted) if j == i[1]]:
			plt.barh(bar_height*i[0] + v_offset + sig_min, width=window, left=sample, height=bar_height, 
				color=labels_style[i[0]], alpha=0.5, linewidth=0)
		"""
	# set the tick markers of both, ground truth and predicted values
	ticks = [-(bar_height*i - bar_height/2) - v_offset + sig_max for i in range(len(labels))]
	ticks.extend([(bar_height*i + bar_height/2) + v_offset + sig_min for i in range(len(labels))])
	plt.yticks(ticks, list(labels.values())*2)

	plt.title(title)

	# add a boxplot
	#plt.subplot(122)
	plt.subplot2grid((1,10), (0,9), colspan=1)
	plt.boxplot(signal)
	plt.title('Box Plot')
	plt.xticks([])
	if ylim:
		plt.ylim(ylim[0], ylim[1])

	plt.subplots_adjust(right=0.98, bottom=0.15, wspace=0.8)
	
	if display == 'plot':
		plt.show()
	elif display == 'save':
		plt.savefig(directory+filename, format='png')
		plt.close()








# TODO allow to set other style parameters other than the color
def plot_conf_mat(conf_mat, labels=[], pred_style=[''], grd_style=[''], grd_text='', pred_text='', 
				mat_title='', plot_title='', ylabel_plot='Percentage of Signal', plot_xlabel_xy=(1, -10),
				mat_xlabel_xy=(0, 0.6), mat_ylabel_xy=(-0.15, 0.6), bar_width=0.4, bar_offset=0.3, 
				wspace=0.2, hspace=0.2, left=0.1, right=0.95, top=0.95, bottom=0.2, title='', 
				display='plot', filename='', directory='', figsize=()):
	"""
####################################################################################
Plot the Confusion Matrix (Contigency Matrix)
####################################################################################
conf_mat		List of list. The confusion matrix as returned by the function 
				sklearn.metrics.confusion_matrix.
labels			List of string. The list with the labels regarding the same order as
				in the confusion matrix. 
pred_style		List. The color of the predicted values labels
grd_style		List. The color of the ground truth labels
grd_text		String. The text for the ground truth
pred_text		String. The tex for the predicted
mat_title		String. The title of the confusion matrix
plot_title		String. The title of the plot (the right one) with the percentages
ylabel_plot		String. The y label of the plot (the right one)
plot_xlabel_xy	Tuple. The position (x,y) of the xlabel on the plot (the right one)
mat_xlabel_xy	Tuple. The position (x,y) of the xlabel on the confusion matrix
mat_ylabel_xy	Tuple. The position (x,y) of the ylabel on the confusion matrix
bar_width		Number. The width of the bars
bar_offset		Number. The offset of the bars
wpsace			Number. The amount of width reserved for blank space between subplots (plot and matrix). See pyplot.subplots_adjust
hspace			Number. The amount of height reserved for white space between subplots.  See pyplot.subplots_adjust
left			Number. The left side of the subplots of the figure. See pyplot.subplots_adjust
right			Number. The right side of the subplots of the figure. See pyplot.subplots_adjust
top				Number. The top of the subplots of the figure. See pyplot.subplots_adjust
bottom			Number. The bottom of the subplots of the figure. See pyplot.subplots_adjust
tite			String. The plot's title
display			String. 'plot' to display the generated plots, 'subplot' for a subplot graphic, something else to save
				to a file
filename		String. The filename with extension
directory		String. The directory where to save the image, remember to end with a slash '/'
figsize			Tuple. The figure size (widht, height) in inches
####################################################################################

Usage example:

from sklearn.metrics import confusion_matrix
from ts2cn.plot import graphics

lb = 'ab'; signal = data[0]; predicted = oo['events'][0][0]; grd_truth = [lb]*len(predicted); labels_dic = {'ab': 'Abnormal', 'no': 'Normal', 'sp':'Spike', 'sz':'Seizure'}; labels_style = ['red', 'blue', 'green', 'yellow']; imp.reload(graphics); graphics.plot_sig_pred_grd(signal, 'Samples', 'Signal', 'c-', grd_truth, 'Ground Truth',  predicted, 'Predicted', labels_dic, labels_style, labels_pos=-600, window=w_size, bar_height=10, v_offset=0, display='subplot', figsize=(15,5))

conf_mat = confusion_matrix(grd_truth, predicted, labels=['no', 'sz', 'sp', 'ab'])


imp.reload(graphincs); graphics.plot_conf_mat(conf_mat, labels=['no', 'sz', 'sp', 'ab'], labels_style=['red', 'green', 'blue', 'yellow'], grd_text='GRD_TEXT', pred_text='PRED_TEXT', predicted=[], pred_ylabel='', labels={}, labels_style=[], labels_pos=-600, window=1, bar_height=10, v_offset=1, title='Title', display='plot', filename='', directory='', figsize=())

imp.reload(graphics); graphics.plot_conf_mat(conf_mat, labels=['no', 'sz', 'sp', 'ab'], pred_style=['red', 'green', 'blue', 'yellow'], grd_style=['red', 'green', 'blue', 'yellow'], grd_text='GRD_TEXT', pred_text='PRED_TEXT', bar_width=0.2, bar_offset=0.4, title='Title', display='plot', filename='', directory='', figsize=(), wspace=0.2, hspace=0.2, top=0.95, bottom=0.2, right=0.95, left=0.1, mat_title='MAT_TITLE', plot_title='PLOT_TITLE', ylabel_plot='YLABEL_PLOT')
	"""
	import numpy as np	
	imp.reload(plt)
	if figsize:
		plt.figure(figsize=figsize)
	
	data = np.ndarray(shape=(len(conf_mat), len(conf_mat[0])), buffer=np.array(conf_mat), dtype=int)
	columns = labels
	rows = labels

	values = np.arange(0, 105, 10)
	value_increment = 1

	n_rows = len(data)

	index = np.arange(len(columns)) + bar_offset

	# Initialize the vertical-offset for the stacked bar chart.
	y_offset = np.array([0.0] * len(columns))
	
	plt.subplot(121)
	plt.grid(b='on')
	data_perc = data/sum(data)*100
	# Plot bars 
	for row in range(n_rows):
		plt.bar(index, data_perc[row], width=bar_width, bottom=y_offset, color=grd_style[row])
		y_offset = y_offset + data_perc[row]
	
	
	# Add colorfull labels to the plot
	the_table = plt.table(cellText=[columns], cellColours=[pred_style], cellLoc='center', loc='bottom')
	

	# Adjust layout to make room for the table:
	plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom, hspace=hspace, wspace=wspace)

	plt.ylabel(ylabel_plot)
	plt.ylim(values[0], values[-1])
	plt.yticks(values * value_increment, ['%d' % val for val in values])
	plt.minorticks_on()
	plt.xticks(index + bar_width/2, [])
	plt.text(s=pred_text, x=plot_xlabel_xy[0], y=plot_xlabel_xy[1])
	plt.title(plot_title)


	plt.subplot(122)
	
	# Add a table at the bottom of the axes
	the_table = plt.table(cellText=data,
						  rowLabels=rows, #rowLoc='right',
						  rowColours=grd_style,
						  colLabels=columns,
						  colColours=pred_style,
						  loc='center')
	plt.text(s=pred_text, x=mat_xlabel_xy[0], y=mat_xlabel_xy[1], axes=plt.gca())
	plt.text(s=grd_text, x=mat_ylabel_xy[0], y=mat_ylabel_xy[1], rotation=90, axes=plt.gca())
	plt.axis('off')
	plt.title(mat_title)

	plt.suptitle(title)
	
	if display == 'plot':
		plt.show()
	elif display == 'save':
		plt.savefig(directory+filename, format='png')
		plt.close()







def plot_roc_curve(predicted, ground_truth, labels, compl=None, titles=[], aggregate=False, x_label='False Positive Rate', 
				y_label='True Positive Rate', suptitle='', leg_fontsize='small', left=0.125, right=0.9, top=0.9, 
				bottom=0.1, wspace=0.2, hspace=0.2, err_str='Not possible to calculate ROC', display='plot', 
				filename='', directory='', figsize=()):
	"""
####################################################################################
Plot signal along side the predicted values and groundtruth sample by sample
####################################################################################
predicted		List. The list of the predicted values
ground_truth	List. The list of the ground truth values
labels			List. The list of the labels
compl			Dict of boolean. The key is the same elements as the list labels, and the values are boolean,
				True to perform the complement comparison in the respective label.
				e.g., {'no':True, 'sz':False} perform the complement comparison for the first label ('no') but 
				not for the second ('sz'). This would compare the first label with "sample != lb" instead of 
				"sample == lb". 
titles			List of string. The list of the string that represent each label, e.g., 'no' means 'Normal'
aggregate		Boolean. True to aggregate the ROC point of several not probabilistic classifiers, instead of 
				interpolaring the ROC curve, using the predicted score.
x_label			String. The string to display at the x axis
y_label			String. The string to display at the y axis
wpsace			Number. The amount of width reserved for blank space between subplots (plot and matrix). See pyplot.subplots_adjust
hspace			Number. The amount of height reserved for white space between subplots.  See pyplot.subplots_adjust
left			Number. The left side of the subplots of the figure. See pyplot.subplots_adjust
right			Number. The right side of the subplots of the figure. See pyplot.subplots_adjust
top				Number. The top of the subplots of the figure. See pyplot.subplots_adjust
bottom			Number. The bottom of the subplots of the figure. See pyplot.subplots_adjust
suptite			String. The plot's title
leg_fontsize	String. The fontsize of the legend (int or float or {‘xx-small’, ‘x-small’, 
				‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’})
err_str			String. The string to display when it is not possible to calculate the ROC curve
display			String. 'plot' to display the generated plots, 'subplot' for a subplot graphic, something else to save
				to a file
filename		String. The filename with extension
directory		String. The directory where to save the image, remember to end with a slash '/'
figsize			Tuple. The figure size (widht, height) in inches
####################################################################################

Usage example:
from ts2cn.plot import graphics

pred = oo['events'][9]; grd = labels; lbs = ['no', 'sz', 'sp', 'ab']; titles = ['Normal', 'Seizure', 'Spike', 'Abnormal']; imp.reload(graphics); ll=[graphics.plot_roc_curve(predicted=oo['events'][i], ground_truth=grd, labels=lbs, titles=titles, figsize=(), hspace=0.4, leg_fontsize='small', suptitle="Roc Curvers TF=%s" % TFs[i], display='save', filename="roc_TF_%s.png" % TFs[i], directory='/tmp/') for i in range(len(TFs))]

data=tss[0][4]['F'][:8]; imp.reload(mt); pp=datetime.now(); levels=3; w_size=int(40/2**levels); labels = [['sp']*int(len(data[0])/w_size)] * len(data); wavelet='db4'; print(pp); oo = mt.automated_categorization(np.hstack(baseline), data, labels, w_size, wavelet, levels, int(21), spike_amp, clean_event_window, TFs[:4], silent=False); perc_grd = np.array([[sum(sample == lb)/len(sample) for lb in ['no', 'sz', 'sp', 'ab']] for sample in np.array(labels)]); perc_pred = np.array([[sum(sample == lb)/len(sample) for lb in ['no', 'sz', 'sp', 'ab']] for sample in np.array(oo['events'][0])]); perc_pred; print(oo['metrics']);


# aggregate without complement
pred = oo['events'][9]; grd = labels; lbs = ['no', 'sz', 'sp', 'ab']; titles = ['Normal', 'Seizure', 'Spike', 'Abnormal']; imp.reload(graphics); pp=datetime.now(); print(pp); ll=graphics.plot_roc_curve(predicted=oo['events'], ground_truth=[grd]*len(oo['events']), labels=lbs, aggregate=True, compl={'no':False, 'sz':False, 'sp':False, 'ab':False}, titles=titles, figsize=(), hspace=0.4, wspace=0.3, leg_fontsize='small', suptitle="Roc Curvers", display='save', filename="roc_TF.png", directory='/tmp/z/'); print(datetime.now()); print( datetime.now() - pp)

# aggregate with complement
pred = oo['events'][9]; grd = labels; lbs = ['no', 'sz', 'sp', 'ab']; titles = ['All Events', 'Seizure', 'Spike', 'Abnormal']; imp.reload(graphics); pp=datetime.now(); print(pp); ll=graphics.plot_roc_curve(predicted=oo['events'], ground_truth=[grd]*len(oo['events']), labels=lbs, aggregate=True, compl={'no':True, 'sz':False, 'sp':False, 'ab':False}, titles=titles, figsize=(), hspace=0.4, wspace=0.3, leg_fontsize='small', suptitle="Roc Curvers", display='save', filename="roc_TF2.png", directory='/tmp/z/'); print(datetime.now()); print( datetime.now() - pp)

# without aggregate without complement
pred = oo['events'][9]; grd = labels; lbs = ['no', 'sz', 'sp', 'ab']; titles = ['Normal', 'Seizure', 'Spike', 'Abnormal']; imp.reload(graphics); pp=datetime.now(); print(pp); ll=[graphics.plot_roc_curve(predicted=oo['events'][i], ground_truth=grd, labels=lbs, compl={'no':False, 'sz':False, 'sp':False, 'ab':False}, titles=titles, figsize=(), hspace=0.4, wspace=0.3, leg_fontsize='small', suptitle="Roc Curvers TF=%s" % TFs[i], display='save', filename="roc_TF_%s.png" % TFs[i], directory='/tmp/t/') for i in range(len(TFs))]; print( datetime.now() - pp)

# without aggregate with complement
pred = oo['events'][9]; grd = labels; lbs = ['no', 'sz', 'sp', 'ab']; titles = ['All Events', 'Seizure', 'Spike', 'Abnormal']; imp.reload(graphics); pp=datetime.now(); print(pp); ll=[graphics.plot_roc_curve(predicted=oo['events'][i], ground_truth=grd, labels=lbs, compl={'no':True, 'sz':False, 'sp':False, 'ab':False}, titles=titles, figsize=(), hspace=0.4, wspace=0.3, leg_fontsize='small', suptitle="Roc Curvers TF=%s" % TFs[i], display='save', filename="roc_TF_%s.png" % TFs[i], directory='/tmp/y/') for i in range(len(TFs))]; print( datetime.now() - pp)


# calculating all confusion_matrixes

from sklearn.metrics import confusion_matrix, auc
lbs = ['no', 'sz', 'sp', 'ab']; mat = np.array([confusion_matrix(y_true=np.hstack(labels), y_pred=np.hstack(oo['events'][i]), labels=lbs) for i in range(len(TFs))])
positive = lbs.index('sz');
tpr = [m[positive][positive] / sum(m[:][positive]) for m in mat]
negative = np.array(range(len(mat[0]))) != positive
fpr = [1- ( sum(sum(m[negative,:][:,negative])) / (sum(sum(m[negative,:][:,negative])) + sum(m[negative,:][:,positive])) ) for m in mat]
plt.plot(fpr, tpr, 'r-*', label="AUC %f.2" % auc(fpr, tpr)); plt.legend(loc='best'); plt.show()
	"""
	
	from sklearn.metrics import roc_curve, auc, confusion_matrix
	import numpy as np
	import imp
	imp.reload(plt)

	if figsize:
		fig = plt.figure(figsize=figsize)
	else:
		fig = plt.figure()
	
	# calculate the number of rows in the plot
	rows = int(np.ceil(np.sqrt(len(labels))))
	columns = int(np.floor(np.sqrt(len(labels))))
	id_plot = 1
	
	# set the default value for compl, the dictionary
	# that specifies which label the comparison should 
	# be complemented
	if compl == None:
		compl = dict([(i, False) for i in labels])
	
	if aggregate:
		# calculate the confusion matrix for each classifier found with given parameters
		# Used when aggregate all classifiers to a single roc curve
		conf_mat = np.array([confusion_matrix(y_true=np.hstack(ground_truth[i]), y_pred=np.hstack(predicted[i]), labels=labels) 
							for i in range(len(predicted))])
	else:
		# calculate the percentage of each label for each instance (sample)
		perc_pred = np.array([[sum(sample != lb)/len(sample) if compl[lb] else sum(sample == lb)/len(sample)
							for lb in labels] for sample in np.array(predicted)])
		perc_grd = np.array([[sum(sample != lb)/len(sample) if compl[lb] else sum(sample == lb)/len(sample) 
							for lb in labels] for sample in np.array(ground_truth)])
	
#	if not aggregate:
#		print('perc_grd'+ str(perc_grd))
#		print('perc_pred'+ str(perc_pred))

	for lb in range(len(labels)):
		# if aggregate, find the point in the ROC space for each classifier,
		# each classifier is regarded as a point and not a whole curve anymore
		if aggregate:
			positive = lb
			negative = np.array(range(len(conf_mat[0]))) != positive

			# if complement a given label, change negative and positive
			if compl[ labels[lb] ]:
				positive = negative
				negative = lb
				tpr = [sum(sum(m[positive, :][:, positive])) / sum(sum(m[positive, :])) for m in conf_mat]
				fpr = [1 - (m[negative, negative] / sum(m[negative, :])) for m in conf_mat]
			else:	
				tpr = [m[positive, positive] / sum(m[positive, :]) for m in conf_mat]
				fpr = [1 - (sum(sum(m[negative, :][:, negative])) / sum(sum(m[negative, :])) ) for m in conf_mat]
			# insert a point at the origin (0,0) and at the (1,1)
			tpr.insert(0, 0)
			fpr.insert(0, 0)
			tpr.append(1)
			fpr.append(1)
			
			
			# sort tpr and fpr
			indexes = np.argsort(fpr)
			fpr = np.array(fpr)[indexes]
			tpr = np.array(tpr)[indexes]
			
		else:
#			print('')
#			print('lb='+ str(lb))
#			print('labels[lb]='+ str(labels[lb]))
#			print(perc_grd[:,lb])
#			print(perc_pred[:,lb])
			fpr, tpr, _ = roc_curve(perc_grd[:,lb], perc_pred[:,lb], pos_label=1)
		
		# check if there is any "nan" in fpr or tpr
		label_str = err_str
		if not np.isnan((fpr, tpr)).any():
			label_str = 'AUC = %0.2f' % auc(fpr, tpr)
		

		plt.subplot(rows, columns, id_plot)
		plt.plot(fpr, tpr, label=label_str)
		plt.plot([0, 1], [0, 1], 'k--')
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.05])
		plt.xlabel(x_label)
		plt.ylabel(y_label)
		plt.title(titles[lb])
		plt.legend(loc="best", fontsize=leg_fontsize)
		
		id_plot += 1
	
	
	plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom, hspace=hspace, wspace=wspace)
	plt.suptitle(suptitle)
	
	if display == 'plot':
		plt.show()
	elif display == 'save':
		plt.savefig(directory+filename, format='png')
		plt.close()
	
	return (fpr, tpr)
	



def recurrence_plot(matrix, display, ths=None, interpolation='none', directory='/tmp/', filename='recurrence.png', 
				xlabel='', ylabel='', title='', figsize=()):
	"""
####################################################################################
 Recurrence Plot
####################################################################################
matrix			2-dimensional Array. 
ths				Number. A threshold to turn the matrix into binary data. Above ths is True, 1, and below is False, 0.
				If None is passed then no threshold is applied. The default is None.
interpolation	String. The type of interpolation. See matplotlib imshow documentation for details
tite			String. The plot's title
display			String. 'plot' to display, 'save' to save to a file
filename		String. The filename with extension
directory		String. The directory where to save the image, remember to end with a slash '/'
figsize			Tuple. The figure size (widht, height) in inches
####################################################################################

usage examples: 

	"""
	m = array(matrix)
	if ths:
		m = m > ths
	if figsize:
		plt.figure(figsize=figsize)
	plt.imshow(m, cmap=plt.cm.gray, interpolation=interpolation)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	
	if display == 'plot':
		plt.show()
	elif display == 'save':
		plt.savefig(directory+filename, format='png')
		plt.close()




def plot_sig_wind_clusters(signal, w_size=1, clusters=[], sig_xlabel='', sig_ylabel='', sig_style='', 
				cls_style='', alpha=0.4, v_offset=0, bar_height=10, title='', display='plot', 
				filename='', directory='', figsize=()):
	"""
####################################################################################
Plot signal along side the predicted values and groundtruth sample by sample
####################################################################################
signal		List. The signal to plot
w_size
clusters
sig_xlabel	String. A label for the signal (xlabel)
sig_ylabel	String. A label for the signal (ylabel)
sig_style	String. The style for the signal
bar_height	Number. The height of the horizontal bars
ylim		Tuple. The ylim min and max
ticks_space	Number. The space between x ticks 
tite		String. The plot's title
display		String. 'plot' to display the generated plots, 'subplot' for a subplot graphic, something else to save
			to a file
filename	String. The filename with extension
directory	String. The directory where to save the image, remember to end with a slash '/'
figsize		Tuple. The figure size (widht, height) in inches
####################################################################################

Usage example:

ts = hstack((ts1[0], ts2[0])); imp.reload(gr); gr.plot_sig_wind_clusters(ts, w_size=256, clusters=r[2][1][1], bar_height=std(ts)*3, cls_style=['black', 'red', 'blue', 'green', 'yellow', 'purple', 'gray', 'orange'], alpha=0.3)

	"""
	imp.reload(plt)

	import numpy as np
	
	if figsize:
		plt.figure(figsize=figsize)
	
	# if no clusters style is passed get and automatic colormap
	if not cls_style:
		cls_style = plt.cm.hot(list(range(len(clusters))))
	
	# plot the horizontal bars
	for i in range(len(clusters)):
		plt.barh(v_offset -bar_height/2, width=w_size, left=i*w_size, height=bar_height, 
				color=cls_style[clusters[i]], alpha=alpha, linewidth=0)
	
	# plot the signal
	plt.plot(signal, sig_style, color='cyan', zorder=-1)
	plt.title(title)
	plt.xlabel(sig_xlabel)
	plt.ylabel(sig_ylabel)
	
	if display == 'plot':
		plt.show()
	elif display == 'save':
		plt.savefig(directory+filename, format='png')
		plt.close()






def plot_windows_stationarity(tss, ts_size, ts_id, w_sizes, cell_test, cell_labels, 
					cell_height=1, cell_width=1, cell_str_color={}, cell_hit_color='lightblue', 
					title='', xlabel='', ylabel='', directory='', filename='', figsize=()):
	"""
####################################################################################

####################################################################################
tss				List. List of dictionaries, each element of the list has a dictionary with the size and an id for
				the time series
w_sizes			List. List with the window sizes used
cell_labels		List. List of lists where each element has a list of the labels for a given window size 
cell_height 	Number. 
cell_width		Number. 
cell_str_color	Dict. Dictionary with the label as key and its value is the color to fill the cell
cell_hit_color	String. The color used when the cell is marked as passed the statistical test
title

filename		String. The filename with extension
directory		String. The directory where to save the image, remember to end with a slash '/'
figsize			Tuple. The figure size (widht, height) in inches
####################################################################################

Usage example:

from matplotlib import pyplot as plt
import imp
from numpy import hstack, array
from datasets import bonn as dt_bonn
from stationarity import bonn

tss = dt_bonn.load()

size = 3; data = hstack([tss['O'], tss['S']])[:size]; labels = [hstack([['no']*len(tss['O'][0]), ['sz']*len(tss['S'][0])])] *size; imp.reload(bonn); hh = bonn.chi_square_analysis('type', 'ID', data, labels, [[]], [[]], [1024], 20, [0.05], './results/tmp/teste/' )

cell_test = [i[2] for i in hh[0][0]['perf_sel'][0]]


from ts2cn.plot import graphics
i=0; imp.reload(graphics); graphics.plot_windows_stationarity([], len(data[i]), 'ID', [1024], cell_test[i], labels[i], 1, 1, {'no': 'black', 'sz':'black'}, 'lightblue', '', '', '', '', '', (10,6))


cell_test2 = [[i[2] for i in b[0]['perf_sel'][0]] for b in hh]
i=0; imp.reload(graphics); graphics.plot_windows_stationarity([], len(data[i]), 'ID', [1024, 2048, 4096], cell_test2, labels[i], 1, 1, {'no': 'black', 'sz':'black'}, 'lightblue', '', '', '', '', '', (10,6))

	"""
	imp.reload(plt)
	#stt = ['AA', 'BB']; bottom=[0,1]; width=[2, 10]; height=[1,2]; left=[0,3]; imp.reload(plt); bars = plt.barh(bottom, width, height, left, color='lightblue'); [plt.annotate(stt[i], (j.xy[0] + width[i]/2, j.xy[1] + height[i]/2), color='blue', fontsize=18, style='normal', weight='bold') for i,j in enumerate(bars)]; plt.show()
	
	default_color = 'white'
	fontsize = 18
	style = 'normal'
	weight = 'bold'
	
	yticks = []
	plt.figure(figsize=figsize)
	for w_i, w_size in enumerate(w_sizes):
		left = [w_size*i for i,j in enumerate(cell_test)]
		color = [cell_hit_color if passed else default_color for i, passed in enumerate(cell_test)]
		width = [w_size * cell_width] * len(cell_test)
		bottom = [w_i*cell_height] * len(cell_test)
		print(('left', left))
		print(('color', color))
		print(('width', width))
		bars = plt.barh(bottom=bottom, width=width, height=cell_height, left=left, color=color)
		for b_i, bar in enumerate(bars):
			plt.annotate(cell_labels[b_i], (bar.xy[0] + width[0]/2, bar.xy[1] + cell_height/2), color=cell_str_color[cell_labels[b_i]], fontsize=fontsize, style=style, weight=weight)
		
		yticks.append("%s - %d"%(ts_id, w_size))
	
	
	plt.yticks(arange(0.5, len(yticks) +1), yticks)
	plt.xticks(range(ts_size), range(ts_size))
	plt.show()
	plt.close()
	
	
