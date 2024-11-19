#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import scipy
import h5py as h5
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

phalanges = np.array([2, 4, 8, 16, 32])
phalanges_labels = ['I-3', 'I-2', 'I-1', 'T-1', 'I-2*']

# Paths
for i_sub in range(13):
	path = '/home/chrpfe/Documents/21099_opm/Phalanges/sub-'+str(i_sub).zfill(2)
	opm_file = 'sub'+str(i_sub).zfill(2)+'_opm_ica_ds.mat'
	meg_file = 'sub'+str(i_sub).zfill(2)+'_meg_ica_ds.mat'

	# --- Read OPM data ---------------------------------
	# OPM data
	f=h5.File(path+opm_file)

	# Read trials
	tmp = f['/opm_ica_ds/trial']
	n_trials = tmp.size
	(n_samples,n_channels) = f.get(tmp[1,0]).shape
	X = np.zeros((n_trials,n_channels,n_samples))
	for i in range(n_trials):
	    X[i,:,:] = np.transpose(np.array(f.get(tmp[i,0])))

	# Read time
	tmp = f['/opm_ica_ds/time']
	n_trials = tmp.size
	n_samples = f.get(tmp[1,0]).size
	time = np.zeros((n_samples,))
	time[:,] = np.transpose(np.array(f.get(tmp[i,0])))

	# Read labels
	tmp = f['/opm_ica_ds/label']
	channel_labels = np.zeros((tmp.size,), dtype='<U7')
	for i in range(tmp.size):
	    channel_labels[i] = np.array(bytes(f.get(tmp[0,i])[:]).decode('utf-16'))

	# Find channels of interest  
	opm_channels = np.zeros((n_channels,), dtype='bool')
	opmeeg_channels = np.zeros((n_channels,), dtype='bool')
	for i in range(channel_labels.size):
	    if 'bz' in channel_labels[i]:
	        opm_channels[i] = 1;
	    if 'EEG' in channel_labels[i]:
	        opmeeg_channels[i] = 1;

	# MEG data labels
	tmp = f['/opm_ica_ds/trialinfo']
	y_opm = np.array(tmp[0])
	# OPM data (mags)
	X_opm_orig = X[:,opm_channels,:]
	X_opm = X_opm_orig.reshape(len(X_opm_orig), -1)
	X_opm = X_opm / np.std(X_opm)
	# EEG data
	X_opmeeg_orig = X[:,opmeeg_channels,:]
	X_opmeeg = X_opmeeg_orig.reshape(len(X_opmeeg_orig), -1)
	X_opmeeg = X_opmeeg / np.std(X_opmeeg)


	# --- Read MEG data -----------------------------
	# MEG data
	f=h5.File(path+meg_file)

	# Read trials
	tmp = f['/meg_ica_ds/trial']
	n_trials = tmp.size
	(n_samples,n_channels) = f.get(tmp[1,0]).shape
	X = np.zeros((n_trials,n_channels,n_samples))
	for i in range(n_trials):
	    X[i,:,:] = np.transpose(np.array(f.get(tmp[i,0])))

	# Read labels
	tmp = f['/meg_ica_ds/label']
	channel_labels = np.zeros((tmp.size,), dtype='<U7')
	for i in range(tmp.size):
	    channel_labels[i] = np.array(bytes(f.get(tmp[0,i])[:]).decode('utf-16'))

	# Find channels of interest  
	meg_channels = np.zeros((n_channels,), dtype='bool')
	megmag_channels = np.zeros((n_channels,), dtype='bool')
	megeeg_channels = np.zeros((n_channels,), dtype='bool')
	for i in range(channel_labels.size):
	    if 'MEG' in channel_labels[i]:
	        meg_channels[i] = 1;
	        if channel_labels[i][-1] == '1':
	            megmag_channels[i] = 1;
	    if 'EEG' in channel_labels[i]:
	        megeeg_channels[i] = 1;

	# MEG data labels
	tmp = f['/meg_ica_ds/trialinfo']
	y_meg = np.array(tmp[0])
	# MEG data (mags and grads)
	X_meg_orig = X[:,meg_channels,:]
	X_meg = X_meg_orig.reshape(len(X_meg_orig), -1)
	X_meg = X_meg / np.std(X_meg)
	# MAG data (mags only)
	X_megmag_orig = X[:,megmag_channels,:]
	X_megmag = X_megmag_orig.reshape(len(X_megmag_orig), -1)
	X_megmag = X_megmag / np.std(X_megmag)
	# EEG data
	X_megeeg_orig = X[:,megeeg_channels,:]
	X_megeeg = X_megeeg_orig.reshape(len(X_megeeg_orig), -1)
	X_megeeg = X_megeeg / np.std(X_megeeg)

	# Find trials for each phalanges and chose number of trials to use
	n_trials = np.min([np.sum(y_meg==2),np.sum(y_meg==4),
		           np.sum(y_meg==8),np.sum(y_meg==16),
		           np.sum(y_meg==32),np.sum(y_opm==2),
		           np.sum(y_opm==4),np.sum(y_opm==8),
		           np.sum(y_opm==16),np.sum(y_opm==32)])

	print("Phalange: OPM trials / MEG trials")
	i_trial_opm = np.zeros((5,n_trials),dtype='int32')
	for i in range(5):
	    i_trial_opm[i,:] = np.arange(y_opm.size)[y_opm==phalanges[i]][0:n_trials]   
	i_trial_meg = np.zeros((5,n_trials),dtype='int32')
	for i in range(5):
	    i_trial_meg[i,:] = np.arange(y_meg.size)[y_meg==phalanges[i]][0:n_trials]
	    print("%s: %d / %d" % (phalanges_labels[i], sum(y_opm==phalanges[i]), sum(y_meg==phalanges[i])))
	print("N_trials: ",n_trials)


	# --- OPM classification ------------------------------------------------ 
	# Bionmial classification
	opm_accuracy = np.zeros((5,5))
	pipe = make_pipeline(StandardScaler(), LogisticRegression())
	for c1 in range(5):
	    for c2 in range((c1+1),5):
	        trl_sel = np.sort(np.concatenate((i_trial_opm[c1,:], i_trial_opm[c2,:])))
	        y_test = y_opm[trl_sel]
	        X_test = X_opm[trl_sel,:]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        opm_accuracy[c1,c2]  =np.mean(scores_full)
	        opm_accuracy[c2,c1]  =np.mean(scores_full)
	        print("Classification score (%s vs %s): %0.3f (std. %0.3f); n_trials = %d" % (phalanges[c1], phalanges[c2], np.mean(scores_full), np.std(scores_full), y_test.size))

	fig, ax = plt.subplots()
	cax = ax.matshow(opm_accuracy)
	for (i, j), z in np.ndenumerate(opm_accuracy):
	    if i!=j:
	        ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
	fig.colorbar(cax)
	plt.title('OPM: Pairwise classification accuracy')
	xaxis = np.arange(len(phalanges));
	ax.set_xticks(xaxis)
	ax.set_yticks(xaxis)
	ax.set_xticklabels(phalanges_labels)
	ax.set_yticklabels(phalanges_labels)
	plt.savefig(fig_path+'LogReg_OPM_pairwise.png', bbox_inches='tight')

	# Multinomial classification
	trl_all = np.sort(np.concatenate((i_trial_opm[0,:], i_trial_opm[1,:], i_trial_opm[2,:], 
		                          i_trial_opm[3,:], i_trial_opm[4,:])))
	X_train, X_test, y_train, y_test = train_test_split(X_opm[trl_all,:], y_opm[trl_all], random_state=0)
	clf = LogisticRegression()
	pipe.fit(X_train, y_train)
	y_pred = pipe.predict(X_test)
	opm_CM = confusion_matrix(y_test, y_pred, labels=pipe.classes_, normalize='true')
	ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', display_labels=phalanges_labels)
	plt.title('OPM: Multinomial classification confusion matrix')
	plt.savefig(fig_path+'LogReg_OPM_multinomial.png', bbox_inches='tight')

	# Classification over time (one vs. rest)
	opm_accuracy_over_time = np.zeros((5,X_opm_orig.shape[2]))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for ph in range(5):
	    tmp = np.empty((1,), dtype='int32')
	    for i in range(5):
	        if i == ph:
	            tmp = np.concatenate((tmp,i_trial_opm[i,0:(int(np.floor(289/4)*4))]))
	        else:
	            tmp = np.concatenate((tmp,i_trial_opm[i,0:int(np.floor(289/4))]))
	    i_tmp = np.sort(tmp[1:])
	    y_test = y_opm[i_tmp]
	    y_test[y_test!=phalanges[ph]] = 0
	    for t in range(X_opm_orig.shape[2]):
	        #X_test = np.vstack((X_false[:,:,t], X_true[:,:,t]))
	        X_test = X_opm_orig[i_tmp,:,t]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        opm_accuracy_over_time[ph,t]  =np.mean(scores_full)

	fig, ax = plt.subplots()
	for i in range(5):
	    ax.plot(time,np.transpose(opm_accuracy_over_time[i,:]),label=phalanges_labels[i])
	plt.title('OPM: Classification (one-vs-rest) accuracy over time')
	plt.legend()
	plt.ylabel('Classification accuracy')
	plt.xlabel('Time [msec]')
	plt.savefig(fig_path+'LogReg_OPM_OneVsRest_overTime.png', bbox_inches='tight')


	# --- OPM-EEG classification -----------------------------------------
	# Binomial classification
	opmeeg_accuracy = np.zeros((5,5))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for c1 in range(5):
	    for c2 in range((c1+1),5):
	        trl_sel = np.sort(np.concatenate((i_trial_opm[c1,:], i_trial_opm[c2,:])))
	        y_test = y_opm[trl_sel]
	        X_test = X_opmeeg[trl_sel,:]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        opmeeg_accuracy[c1,c2]  =np.mean(scores_full)
	        opmeeg_accuracy[c2,c1]  =np.mean(scores_full)
	        print("Classification score (%s vs %s): %0.3f (std. %0.3f)" % (phalanges[c1], phalanges[c2], np.mean(scores_full), np.std(scores_full)))

	fig, ax = plt.subplots()
	cax = ax.matshow(opmeeg_accuracy)
	for (i, j), z in np.ndenumerate(opmeeg_accuracy):
	    if i!=j:
	        ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
	fig.colorbar(cax)
	plt.title('OPM-EEG: Pairwise classification accuracy')
	xaxis = np.arange(len(phalanges));
	ax.set_xticks(xaxis)
	ax.set_yticks(xaxis)
	ax.set_xticklabels(phalanges_labels)
	ax.set_yticklabels(phalanges_labels)
	plt.savefig(fig_path+'LogReg_OPMEEG_pairwise.png', bbox_inches='tight')

	# Multinomial classification
	trl_all = np.sort(np.concatenate((i_trial_opm[0,:], i_trial_opm[1,:], i_trial_opm[2,:], 
		                          i_trial_opm[3,:], i_trial_opm[4,:])))
	X_train, X_test, y_train, y_test = train_test_split(X_opmeeg[trl_all,:], y_opm[trl_all], random_state=0)
	pipe.fit(X_train, y_train)
	y_pred = pipe.predict(X_test)
	opmeeg_CM = confusion_matrix(y_test, y_pred, labels=pipe.classes_, normalize='true')
	ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', display_labels=phalanges_labels)
	plt.title('OPM-EEG: Multinomial classification confusion matrix')
	plt.savefig(fig_path+'LogReg_OPMEEG_multinomial.png', bbox_inches='tight')

	# Classification over time (one vs. rest)
	opmeeg_accuracy_over_time = np.zeros((5,X_opmeeg_orig.shape[2]))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for ph in range(5): 
	    tmp = np.empty((1,), dtype='int32')
	    for i in range(5):
	        if i == ph:
	            tmp = np.concatenate((tmp,i_trial_opm[i,0:(int(np.floor(289/4)*4))]))
	        else:
	            tmp = np.concatenate((tmp,i_trial_opm[i,0:int(np.floor(289/4))]))
	    i_tmp = np.sort(tmp[1:])
	    y_test = y_opm[i_tmp]
	    y_test[y_test!=phalanges[ph]] = 0
	    for t in range(X_opmeeg_orig.shape[2]):
	        #X_test = np.vstack((X_false[:,:,t], X_true[:,:,t]))
	        X_test = X_opmeeg_orig[i_tmp,:,t]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        opmeeg_accuracy_over_time[ph,t]  =np.mean(scores_full)

	fig, ax = plt.subplots()
	for i in range(5):
	    ax.plot(time,np.transpose(opmeeg_accuracy_over_time[i,:]),label=phalanges_labels[i])
	plt.title('OPM-EEG: Classification (one-vs-rest) accuracy over time')
	plt.legend()
	plt.ylabel('Classification accuracy')
	plt.xlabel('Time [msec]')
	plt.savefig(fig_path+'LogReg_OPMEEG_OneVsRest_overTime.png', bbox_inches='tight')


	# --- MEG classification -------------------------------------------------
	# Binomial classification
	meg_accuracy = np.zeros((5,5))
	pipe = make_pipeline(StandardScaler(), LogisticRegression())
	for c1 in range(5):
	    for c2 in range((c1+1),5):
	        trl_sel = np.sort(np.concatenate((i_trial_meg[c1,:], i_trial_meg[c2,:])))
	        y_test = y_meg[trl_sel]
	        X_test = X_meg[trl_sel,:]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        meg_accuracy[c1,c2]  =np.mean(scores_full)
	        meg_accuracy[c2,c1]  =np.mean(scores_full)
	        print("Classification score (%s vs %s): %0.3f (std. %0.3f)" % (phalanges[c1], phalanges[c2], np.mean(scores_full), np.std(scores_full)))

	fig, ax = plt.subplots()
	cax = ax.matshow(meg_accuracy)
	for (i, j), z in np.ndenumerate(meg_accuracy):
	    if i!=j:
	        ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
	fig.colorbar(cax)
	plt.title('MEG: Pairwise classification accuracy')
	xaxis = np.arange(len(phalanges));
	ax.set_xticks(xaxis)
	ax.set_yticks(xaxis)
	ax.set_xticklabels(phalanges_labels)
	ax.set_yticklabels(phalanges_labels)
	plt.savefig(fig_path+'LogReg_MEG_pairwise.png', bbox_inches='tight')

	# Multinomial classification
	trl_all = np.sort(np.concatenate((i_trial_meg[0,:], i_trial_meg[1,:], i_trial_meg[2,:], 
		                          i_trial_meg[3,:], i_trial_meg[4,:])))
	X_train, X_test, y_train, y_test = train_test_split(X_meg[trl_all,:], y_meg[trl_all], random_state=0)
	pipe.fit(X_train, y_train)
	y_pred = pipe.predict(X_test)
	meg_CM = confusion_matrix(y_test, y_pred, labels=pipe.classes_, normalize='true')
	ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', display_labels=phalanges_labels)
	plt.title('MEG: Multinomial classification confusion matrix')
	plt.savefig(fig_path+'LogReg_MEG_multinomial.png', bbox_inches='tight')

	# Classification over time (one vs. rest)
	meg_accuracy_over_time = np.zeros((5,X_meg_orig.shape[2]))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for ph in range(5): 
	    tmp = np.empty((1,), dtype='int32')
	    for i in range(5):
	        if i == ph:
	            tmp = np.concatenate((tmp,i_trial_meg[i,0:(int(np.floor(289/4)*4))]))
	        else:
	            tmp = np.concatenate((tmp,i_trial_meg[i,0:int(np.floor(289/4))]))
	    i_tmp = np.sort(tmp[1:])
	    y_test = y_meg[i_tmp]
	    y_test[y_test!=phalanges[ph]] = 0
	    for t in range(X_meg_orig.shape[2]):
	        #X_test = np.vstack((X_false[:,:,t], X_true[:,:,t]))
	        X_test = X_meg_orig[i_tmp,:,t]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        meg_accuracy_over_time[ph,t]  =np.mean(scores_full)

	fig, ax = plt.subplots()
	for i in range(5):
	    ax.plot(time,np.transpose(meg_accuracy_over_time[i,:]),label=phalanges_labels[i])
	plt.title('MEG: Classification (one-vs-rest) accuracy over time')
	plt.legend()
	plt.ylabel('Classification accuracy')
	plt.xlabel('Time [msec]')
	plt.savefig(fig_path+'LogReg_MEG_OneVsRest_overTime.png', bbox_inches='tight')


	# --- MEG-MAG classification -------------------------------------
	# Binomial classification
	megmag_accuracy = np.zeros((5,5))
	pipe = make_pipeline(StandardScaler(), LogisticRegression())
	for c1 in range(5):
	    for c2 in range((c1+1),5):
	        trl_sel = np.sort(np.concatenate((i_trial_meg[c1,:], i_trial_meg[c2,:])))
	        y_test = y_meg[trl_sel]
	        X_test = X_megmag[trl_sel,:]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        megmag_accuracy[c1,c2]  =np.mean(scores_full)
	        megmag_accuracy[c2,c1]  =np.mean(scores_full)
	        print("Classification score (%s vs %s): %0.3f (std. %0.3f)" % (phalanges[c1], phalanges[c2], np.mean(scores_full), np.std(scores_full)))

	fig, ax = plt.subplots()
	cax = ax.matshow(megmag_accuracy)
	for (i, j), z in np.ndenumerate(megmag_accuracy):
	    if i!=j:
	        ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
	fig.colorbar(cax)
	plt.title('MEG-MAG: Pairwise classification accuracy')
	xaxis = np.arange(len(phalanges));
	ax.set_xticks(xaxis)
	ax.set_yticks(xaxis)
	ax.set_xticklabels(phalanges_labels)
	ax.set_yticklabels(phalanges_labels)
	plt.savefig(fig_path+'LogReg_MEGMAG_pairwise.png', bbox_inches='tight')

	# Multinomial classification
	trl_all = np.sort(np.concatenate((i_trial_meg[0,:], i_trial_meg[1,:], i_trial_meg[2,:], 
		                          i_trial_meg[3,:], i_trial_meg[4,:])))
	X_train, X_test, y_train, y_test = train_test_split(X_megmag[trl_all,:], y_meg[trl_all], random_state=0)
	pipe.fit(X_train, y_train)
	y_pred = pipe.predict(X_test)
	megmag_CM = confusion_matrix(y_test, y_pred, labels=pipe.classes_, normalize='true')
	ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', display_labels=phalanges_labels)
	plt.title('MEG-MAG: Multinomial classification confusion matrix')
	plt.savefig(fig_path+'LogReg_MEGMAG_multinomial.png', bbox_inches='tight')

	# Classification over time (one vs. rest)
	megmag_accuracy_over_time = np.zeros((5,X_megmag_orig.shape[2]))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for ph in range(5): 
	    tmp = np.empty((1,), dtype='int32')
	    for i in range(5):
	        if i == ph:
	            tmp = np.concatenate((tmp,i_trial_meg[i,0:(int(np.floor(289/4)*4))]))
	        else:
	            tmp = np.concatenate((tmp,i_trial_meg[i,0:int(np.floor(289/4))]))
	    i_tmp = np.sort(tmp[1:])
	    y_test = y_meg[i_tmp]
	    y_test[y_test!=phalanges[ph]] = 0
	    for t in range(X_megmag_orig.shape[2]):
	        #X_test = np.vstack((X_false[:,:,t], X_true[:,:,t]))
	        X_test = X_megmag_orig[i_tmp,:,t]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        megmag_accuracy_over_time[ph,t]  =np.mean(scores_full)

	fig, ax = plt.subplots()
	for i in range(5):
	    ax.plot(time,np.transpose(megmag_accuracy_over_time[i,:]),label=phalanges_labels[i])
	plt.title('MEG-MAG: Classification (one-vs-rest) accuracy over time')
	plt.legend()
	plt.ylabel('Classification accuracy')
	plt.xlabel('Time [msec]')
	plt.savefig(fig_path+'LogReg_MEGMAG_OneVsRest_overTime.png', bbox_inches='tight')


	# --- MEG-EEG classification -----------------------------------------
	# Binomial classification
	megeeg_accuracy = np.zeros((5,5))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for c1 in range(5):
	    for c2 in range((c1+1),5):
	        trl_sel = np.sort(np.concatenate((i_trial_meg[c1,:], i_trial_meg[c2,:])))
	        y_test = y_meg[trl_sel]
	        X_test = X_megeeg[trl_sel,:]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        megeeg_accuracy[c1,c2]  =np.mean(scores_full)
	        megeeg_accuracy[c2,c1]  =np.mean(scores_full)
	        print("Classification score (%s vs %s): %0.3f (std. %0.3f)" % (phalanges[c1], phalanges[c2], np.mean(scores_full), np.std(scores_full)))

	fig, ax = plt.subplots()
	cax = ax.matshow(megeeg_accuracy)
	for (i, j), z in np.ndenumerate(megeeg_accuracy):
	    if i!=j:
	        ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
	fig.colorbar(cax)
	plt.title('MEG-EEG: Pairwise classification accuracy')
	xaxis = np.arange(len(phalanges));
	ax.set_xticks(xaxis)
	ax.set_yticks(xaxis)
	ax.set_xticklabels(phalanges_labels)
	ax.set_yticklabels(phalanges_labels)
	plt.savefig(fig_path+'LogReg_MEGEEG_pairwise.png', bbox_inches='tight')

	# Multinomial classification
	trl_all = np.sort(np.concatenate((i_trial_meg[0,:], i_trial_meg[1,:], i_trial_meg[2,:], 
		                          i_trial_meg[3,:], i_trial_meg[4,:])))
	X_train, X_test, y_train, y_test = train_test_split(X_megeeg[trl_all,:], y_meg[trl_all], random_state=0)
	pipe.fit(X_train, y_train)
	y_pred = pipe.predict(X_test)
	megeeg_CM = confusion_matrix(y_test, y_pred, labels=pipe.classes_, normalize='true')
	ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', display_labels=phalanges_labels)
	plt.title('MEG-EEG: Multinomial classification confusion matrix')
	plt.savefig(fig_path+'LogReg_MEGEEG_multinomial.png', bbox_inches='tight')

	# Classification over time (one vs. rest)
	megeeg_accuracy_over_time = np.zeros((5,X_megeeg_orig.shape[2]))
	pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=500))
	for ph in range(5): 
	    tmp = np.empty((1,), dtype='int32')
	    for i in range(5):
	        if i == ph:
	           tmp = np.concatenate((tmp,i_trial_meg[i,0:(int(np.floor(289/4)*4))]))
	        else:
	            tmp = np.concatenate((tmp,i_trial_meg[i,0:int(np.floor(289/4))]))
	    i_tmp = np.sort(tmp[1:])
	    y_test = y_meg[i_tmp]
	    y_test[y_test!=phalanges[ph]] = 0
	    for t in range(X_megeeg_orig.shape[2]):
	        #X_test = np.vstack((X_false[:,:,t], X_true[:,:,t]))
	        X_test = X_megeeg_orig[i_tmp,:,t]
	        scores_full = cross_val_score(pipe, X_test, y_test, cv=10)
	        megeeg_accuracy_over_time[ph,t]  =np.mean(scores_full)

	fig, ax = plt.subplots()
	for i in range(5):
	    ax.plot(time,np.transpose(megeeg_accuracy_over_time[i,:]),label=phalanges_labels[i])
	plt.title('MEG-EEG: Classification (one-vs-rest) accuracy over time')
	plt.legend()
	plt.ylabel('Classification accuracy')
	plt.xlabel('Time [msec]')
	plt.savefig(fig_path+'LogReg_MEGEEG_OneVsRest_overTime.png', bbox_inches='tight')

	# Save scores
	np.savez(path+'pairwise', opm_ac=opm_accuracy, opmeeg_ac=opmeeg_accuracy, meg_ac=meg_accuracy, megmag_ac=megmag_accuracy, megeeg_ac=megeeg_accuracy)
	np.savez(path+'multinomial', opm_cm=opm_CM, opmeeg_cm=opmeeg_CM, meg_cm=meg_CM, megmag_cm=megmag_CM, megeeg_cm=megeeg_CM)
	np.savez(path+'OneVsRest_overTime', opm_acT=opm_accuracy_over_time, opmeeg_acT=opmeeg_accuracy_over_time, meg_acT=meg_accuracy_over_time, megmag_acT=megmag_accuracy_over_time, megeeg_acT=megeeg_accuracy_over_time)



