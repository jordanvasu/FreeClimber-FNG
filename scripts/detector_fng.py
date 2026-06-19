#!/usr/bin/env python
# -*- coding: utf-8 -*- 

## File name : detector.py
## Created by: Adam N. Spierer
## Date      : December 2020
## Purpose   : Script contains main functions used in FreeClimber package, as well as added functionality

version = '0.4.0'
publication = False

import ast
import gc
import os
import sys
import time
import ffmpeg

import numpy as np
import pandas as pd
import trackpy as tp
import trackpy.predict
import subprocess as sp
from scipy.stats import linregress
from scipy.signal import find_peaks,peak_prominences

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

import tortuosity as _tortuosity

## Issue with 'SettingWithCopyWarning' in step_3
pd.options.mode.chained_assignment = None  # default='warn'

class detector(object):
    '''Particle detection platform for identifying the group climbing velocity of a 
    group of flies (or particles) in a Drosophila negative geotaxis (climbing) assay.
    This platform is designed to take videos of varying background homogeneity and 
    applying a background subtraction step before extracting the x,y,time coordinates of
    spots/flies. Climbing velocities are calculated as the most linear portion of a user-defined
    subset of frames by vial (vertical divisions of evenly spaced bins from the min/max
    X-range.
    '''
    def __init__(self, video_file, config_file = None, gui = False, variables = None, debug = False, **kwargs):
        '''Initializing detector object
        ----
        Inputs:
          video_file (str): Path to video file to process
          config_file (str): Path to configuration (.cfg) file
          gui (bool): GUI-specific functions
          variables: None if importing from a file, or list if doing so manually
          debug (bool): Prints out each function as it runs.
          **kwargs: Keyword arguments that are unspecified but can be passed to various plot functions
        ----
        Returns:
          None -- variables are saved to the detector object
        '''
        self.debug = debug
        if self.debug: print('detector.__init__')

        ## Physical vial IDs, populated only from a vials.txt 'id =' line
        ## (see apply_vials_sidecar). None means "label vials by position".
        self.vial_labels = None

        self.config_file = config_file
        self.video_file = self.check_video(video_file)
        
        ## Load variables
        if gui:
            self.load_for_gui(variables = variables)
        elif  not gui:
            self.load_for_main(config_file = config_file)
            ## Per-folder override: a 'vials.txt' beside the video takes
            ## precedence over the cfg's 'vials' (batch convenience for
            ## folders whose videos have a different vial count).
            self.apply_vials_sidecar(video_file)
        self.check_variable_formats()

        ## Setting a color map
        self.vial_color_map = cm.jet
        
        ## Create a conversion factor
        if self.convert_to_cm_sec: self.conversion_factor = self.pixel_to_cm / self.frame_rate
        else: self.conversion_factor = 1

        print('')
        self.specify_paths_details(video_file)
        self.image_stack = self.video_to_array(video_file,loglevel='panic')
        return

    ## Loading functions    
    def load_for_gui(self,variables):
        '''Loads experimental and detections variables to the detector object for GUI
        ----
        Inputs:
          variables (list): List of variables found in the configuration file, order is irrelevant
        ----
        Returns:
          None -- variables are saved to the detector object'''    
        if self.debug: print('detector.load_for_gui')
        self.config = None
        self.path_project = None
        
        ## Exit program if no variables are passed
        if variables == None:
            print('\n\nExiting program. No variable list loaded')
            raise SystemExit
        
        ALLOWED_KEYS = {
            'x', 'y', 'w', 'h', 'check_frame', 'blank_0', 'blank_n',
            'crop_0', 'crop_n', 'threshold', 'diameter', 'minmass',
            'maxsize', 'ecc_low', 'ecc_high', 'vials', 'window',
            'pixel_to_cm', 'frame_rate', 'vial_id_vars', 'outlier_TB',
            'outlier_LR', 'naming_convention', 'path_project',
            'file_suffix', 'convert_to_cm_sec', 'trim_outliers',
            'fng_enabled', 'fng_smooth_window', 'fng_climb_thresh',
            'fng_fall_thresh', 'fng_min_gap', 'fng_recovery_thresh',
            'analysis_mode', 'link_search_range', 'link_memory',
            'link_predictor', 'link_min_track_length',
            'tortuosity_enabled', 'tortuosity_smoothing_window',
            'tortuosity_velocity_threshold', 'tortuosity_bout_min_frames',
            'tortuosity_bout_min_displacement',
        }

        ## Pass imported variables to the detector object
        if self.debug: print('detector.load_for_gui: --------variables--------')
        for item in variables:
            if self.debug: print('detector.load_for_gui:', item)
            if not item.startswith((' ', '\t', '\n')):
                if '=' not in item:
                    continue
                key, _, val_str = item.partition('=')
                key = key.strip()
                if key not in ALLOWED_KEYS:
                    continue
                try:
                    setattr(self, key, ast.literal_eval(val_str.strip()))
                except Exception:
                    ## ast.literal_eval rejects Windows paths (backslash escape
                    ## sequences, e.g. \\U in C:\\Users) and other non-literal
                    ## strings; fall back to the raw value with any surrounding
                    ## quotes stripped rather than dropping the key.
                    setattr(self, key, val_str.strip().strip('"').strip("'"))
                    if self.debug:
                        print('detector.load_for_gui: kept ( %s ) as raw string' % item)
        return
        
    def load_for_main(self, config_file = None):
        '''Loads experimental and detections variables to the detector object for command 
        line interface.
        ----
        Inputs:
          config_file (str): Path to configuration (.cfg) file
        ----
        Returns:
          None -- variables are saved to the detector object''' 
        if self.debug: print('detector.load_for_main')
        
        ## Read in lines from configuration (.cfg) file
        try:
            with open(config_file,'r') as f:
                variables = f.readlines()
            f.close()
            
            ALLOWED_KEYS = {
                'x', 'y', 'w', 'h', 'check_frame', 'blank_0', 'blank_n',
                'crop_0', 'crop_n', 'threshold', 'diameter', 'minmass',
                'maxsize', 'ecc_low', 'ecc_high', 'vials', 'window',
                'pixel_to_cm', 'frame_rate', 'vial_id_vars', 'outlier_TB',
                'outlier_LR', 'naming_convention', 'path_project',
                'file_suffix', 'convert_to_cm_sec', 'trim_outliers',
                'fng_enabled', 'fng_smooth_window', 'fng_climb_thresh',
                'fng_fall_thresh', 'fng_min_gap', 'fng_recovery_thresh',
                'analysis_mode', 'link_search_range', 'link_memory',
                'link_predictor', 'link_min_track_length',
                'tortuosity_enabled', 'tortuosity_smoothing_window',
                'tortuosity_velocity_threshold', 'tortuosity_bout_min_frames',
                'tortuosity_bout_min_displacement',
            }

            ## Filter, format, and import variables to detector object
            if self.debug: print('detector.load_for_main:  --------variables--------')
            variables = [item.rstrip() for item in variables if not item.startswith(('#', ' ', '\t', '\n'))]

            for item in variables:
                if self.debug: print('detector.load_for_main:', item)
                if '=' not in item:
                    continue
                key, _, val_str = item.partition('=')
                key = key.strip()
                if key not in ALLOWED_KEYS:
                    continue
                try:
                    setattr(self, key, ast.literal_eval(val_str.strip()))
                except Exception:
                    print('detector.load_for_main: !! Could not import ( %s )' % item)
            return

        ## Exit program if issue with the configuration file
        except:
            print('\n\nExiting program. Could not read in file.cfg, but path and suffix are good--likely a formatting issue')
            raise SystemExit
        return

    def apply_vials_sidecar(self, video_file):
        '''Override self.vials (and optionally the vial IDs) from a per-folder
        sidecar file, if present.

        Looks for a file named 'vials.txt' in the same folder as the video. This
        lets a single shared .cfg be reused across folders whose videos have
        different vial counts (e.g. 3 or 4 instead of 5): drop a 'vials.txt'
        containing the integer count into the folder and it takes precedence
        over the 'vials' value in the .cfg for every video in that folder.

        Accepted contents (each on its own line, '#'-comment lines ignored):
          * A vial count, as a bare integer ('3') or 'vials=3' / 'vials: 3'.
          * An optional ID line naming which physical vials remain, left to
            right, as comma-separated values: 'id = 33, 34, 35' (also accepts
            'ids', and ':' instead of '='). When present, these IDs label the
            vials in every per-vial output (slopes/fng/tracks/tortuosity) and
            the count is taken from how many IDs are listed -- so an 'id' line
            alone is enough; a separate count line is optional and, if it
            disagrees, the ID count wins.

        A missing or unparseable file leaves self.vials unchanged so the batch
        never breaks.
        ----
        Inputs:
          video_file (str): Video file path
        ----
        Returns:
          None -- updates self.vials (and self.vial_labels) in place when a
                  valid sidecar is found'''
        if self.debug: print('detector.apply_vials_sidecar')

        sidecar = os.path.join(os.path.split(video_file)[0], 'vials.txt')
        if not os.path.isfile(sidecar):
            return

        try:
            with open(sidecar) as f:
                lines = [ln.strip() for ln in f]
            count = None
            labels = None
            for ln in lines:
                if ln == '' or ln.startswith('#'):
                    continue

                ## Identify the key (if any) so 'id'/'ids' lines are routed to
                ## the label parser and everything else is treated as a count.
                key = ''
                value_str = ln
                for sep in ('=', ':'):
                    if sep in ln:
                        key, _, value_str = ln.partition(sep)
                        key = key.strip().lower()
                        value_str = value_str.strip()
                        break

                if key in ('id', 'ids'):
                    labels = self._parse_vial_labels(value_str)
                elif count is None:
                    count = int(value_str)

            if labels:
                ## IDs are authoritative for the count: a video that kept only
                ## vials 33-35 has exactly len(labels) vials regardless of any
                ## (possibly stale) count line.
                if count is not None and count != len(labels):
                    print('!! vials.txt: count (%s) disagrees with %s id(s); using id count :: %s'
                          % (count, len(labels), sidecar))
                self.vial_labels = labels
                print('-- vials.txt override: vials %s -> %s, ids = %s :: %s'
                      % (self.vials, len(labels), labels, sidecar))
                self.vials = len(labels)
            elif count is not None:
                print('-- vials.txt override: vials %s -> %s :: %s' % (self.vials, count, sidecar))
                self.vials = count
            else:
                print('!! vials.txt found but empty, keeping cfg vials = %s :: %s' % (self.vials, sidecar))
        except (ValueError, OSError) as e:
            print('!! Could not parse vials.txt (%s), keeping cfg vials = %s :: %s' % (e, self.vials, sidecar))
        return

    @staticmethod
    def _parse_vial_labels(value_str):
        '''Parse the comma-separated values of a vials.txt 'id =' line into a
        list of vial labels. Numeric entries are kept as ints (so they sort and
        print like vial numbers); any non-numeric entry is kept as a string.
        Returns None when nothing parseable is present.'''
        items = [tok.strip() for tok in value_str.split(',') if tok.strip() != '']
        if not items:
            return None
        labels = []
        for tok in items:
            try:
                labels.append(int(tok))
            except ValueError:
                labels.append(tok)
        return labels

    def _vial_label(self, v):
        '''Translate a positional vial index (1..self.vials) into the physical
        vial ID supplied via a vials.txt 'id =' line, when present. Returns the
        input unchanged when no labels are configured or the index is out of
        range (e.g. the 'all' cohort aggregate).'''
        labels = getattr(self, 'vial_labels', None)
        if labels:
            try:
                iv = int(v)
            except (TypeError, ValueError):
                return v
            if 1 <= iv <= len(labels):
                return labels[iv - 1]
        return v

    def _relabel_vial_col(self, df):
        '''Return a copy of df whose 'vial' column is translated to physical
        vial IDs (see _vial_label), for writing user-facing CSVs. No-op (returns
        df unchanged) when no labels are configured or df has no 'vial' column,
        so the internal positional numbering used for binning/plotting is never
        disturbed.'''
        if not getattr(self, 'vial_labels', None):
            return df
        if 'vial' not in getattr(df, 'columns', []):
            return df
        out = df.copy()
        out['vial'] = out['vial'].map(self._vial_label)
        return out

    ## Checking video path is valid
    def check_video(self,video_file=None):
        '''Checking video path is valid, exiting if not.
        Input:
          video_file (str): Video file path
        ----
        Returns:
          video_file (str): Video file path'''
        if self.debug: print('detector.check_video...', end='')
        
        ## Check if file path is valid, exit if not
        if os.path.isfile(video_file):
            return video_file
        else:
            print('\n\nExiting program. Invalid path to video file: ',video_file)
            raise SystemExit

    ## Specifying variables
    def specify_paths_details(self,video_file):
        '''Takes in the video file and other imported variables and parses them as needed.
        ----
        Inputs:
          video_file (str): Video file path
        ----
        Returns:
          None -- Passes parsed variables back to detector object
        '''
        if self.debug: print('detector.specify_paths_details')
        
        ## Set file and path names
        folder,name = os.path.split(video_file)
        self.name = name[:-5]
        self.name_nosuffix = '.'.join(video_file.split('.')[:-1])
        
        ## Defining final file names and destinations
        file_names = ['data','filtered','diagnostic','slope']
        file_suffixes = ['.raw.csv','.filtered.csv','.diagnostic.png','.slopes.csv']
        for item, jtem in zip(file_names, file_suffixes):
            file_path = self.name_nosuffix + jtem
            if self.debug: print('detector.specify_paths_details: self.path_' + item + "='" + file_path + "'")
            setattr(self, 'path_' + item, file_path)

        ## Project folder specific paths
        if self.path_project == None: self.path_project = os.path.join(folder,self.name + '.cfg')

        ## For future release
#         self.path_review_diagnostic = self.path_project + '/review_at_R_%s/review.log' %self.review_R
        
        ## Extracting file details and naming individual vials
        self.file_details = dict(zip(self.naming_convention.split('_'),self.name.split('_')))
        self.experiment_details = self.name.split('_')
        self.experiment_details[-1] = '.'.join(self.experiment_details[-1].split('.')[:-1])
        self.vial_ID = self.experiment_details[:self.vial_id_vars]
        
        ## Creating a list of colors for plotting
        self.color_list = [self.vial_color_map(i) for i in np.linspace(0,1,self.vials)]
        return

    ## Checking to make sure variables are entered properly...still more to include
    def check_variable_formats(self):
        '''Checks to make sure at least some of the variables input formatted properly'''
        if self.debug: print('detector.check_variable_formats')
        
        ## Vial number
        if self.vials < 1: 
            print('!! Issue with vials: now = 1')
            self.vials = 1
        
        ## Spot diameter must be odd number
        if self.diameter % 2 == 0:
            print('!! Issue with diameter: was %s, now %s' %(self.diameter,self.diameter+1))
            self.diameter += 1
            
        ## Frame rate must be greater than 0
        if self.frame_rate <= 0:
            print('!! Issue with frame_rate: was %s, now 1' %(self.frame_rate))
            self.frame_rate = 1
        
        ## Background blank cannot be greater than the frame 
        if self.blank_0 < self.crop_0:
            self.blank_0 = self.crop_0
        
        if self.blank_n > self.crop_n:
            self.blank_n = self.crop_n

        ## Window size vs. frames to test
        if (self.crop_n - self.crop_0) < self.window:
#             print('!! Issue with window size (%s) being greater than frames (%s). Window size set to 80 percent of desired frames (%s)' %(self.window,self.crop_n - self.crop_0,(.8 * (self.crop_n - self.crop_0))))
            self.window = (self.crop_n - self.crop_0) * 0.8
		
		## blank vs. crop frames
        if self.blank_0 < self.crop_0:
            print('!! Issue with blank frames vs. crop frames. Setting blank_0 (%s) = crop_n (%s)' % (self.blank_0,self.crop_0))
            self.blank_0 = self.crop_0
        if self.blank_n > self.crop_n:
            print('!! Issue with blank frames vs. crop frames. Setting blank_n (%s) = crop_n (%s)' % (self.blank_n,self.crop_n))
            self.blank_n = self.crop_n
        
        ## Check frame is still valid
        if self.check_frame < self.crop_0:
#             print('!! Issue with check_frame < crop_0 (min. cropped frame). Now, check_frame = crop_0 = %s' %self.check_frame)
            self.check_frame = self.crop_0
        if self.check_frame > self.crop_n:
#             print('!! Issue with check_frame > crop_n (max cropped frame). Now, check_frame = crop_n = %s' %self.check_frame)
            self.check_frame = self.crop_n
        return

    ## Video processing functions
    def video_to_array(self, file, **kwargs):
        '''Converts video into an nd-array using ffmpeg-python module.
        ----
        Inputs:
          file (str): Path to video file
          kwargs: Can be used with the ffmpeg.output() argument.
        ----
        Returns:
          image_stack (nd-array): nd-array of the video
        '''
        if self.debug: print('detector.video_to_array')
        
        ## Extracting video meta-data
        try:
            try:
                probe = ffmpeg.probe(file)
            except:
                print('!! Could not read %s into FreeClimber. Likely due to unacceptable video file or FFmpeg not installed' % file)
                raise SystemExit
            video_info = next(x for x in probe['streams'] if x['codec_type'] == 'video')
            self.width = int(video_info['width'])
            self.height = int(video_info['height'])
        except:
            print('!! Could not read in video file metadata')
        
        ## Converting video to nd-array    
        try:
            out,err = (ffmpeg
                       .input(file)
                       .output('pipe:',format='rawvideo', pix_fmt='rgb24',**kwargs)
                       .run(capture_stdout=True))
            self.n_frames = int(len(out)/self.height/self.width/3)
            image_stack = np.frombuffer(out, np.uint8).reshape([-1, self.height, self.width, 3])
        except:
            print('!! Could not read in video file to an array. Error message (if any):', err)

        return image_stack

    def crop_and_grayscale(self,video_array,
                         x = 0 ,x_max = None,
                         y = 0 ,y_max = None,
                         first_frame = None,
                         last_frame = None,
                         grayscale = True):
        '''Crops imported video array to region of interest and converts it to grayscale
        ----
        Inputs:
          video_array (nd-array): image_stack generated from video_to_array function
          x (int): left-most x-position
          x_max (int): right-most x-position
          y (int): lowest y-position
          y_max (int): highest y-position
          first_frame (int): first frame to include
          last_frame (int): last frame to include
          grayscale (bool): True to convert to gray, False to leave in color. 
                                NOTE: Must be in grayscale for FreeClimber, option is 
                                available for functionality beyond FreeClimber
        ----
        Returns:
          clean_stack (nd-array): Cropped and grayscaled (if indicated) video as nd-array'''
        if self.debug: print('detector.crop_and_grayscale')

        ## Conditionals for cropping frames and video length
        if first_frame == None: first_frame = self.crop_0
        if last_frame == None: last_frame = self.crop_n
        if y_max == None: y_max = video_array.shape[2]
        if x_max == None: x_max = video_array.shape[1]
        if self.debug: print('detector.crop_and_grayscale: Cropping from frame %s to %s' % (first_frame,last_frame))
    
        ## Setting only frames and ROI to grayscale
        if grayscale:
            if self.debug: print('detector.crop_and_grayscale: Converting to grayscale & cropping ROI to (%s x %s)' % (x_max-x,y_max-y))
            clean_stack = video_array[first_frame:last_frame,y : y_max,x : x_max,0].astype(np.float64)
            clean_stack *= 0.2989
            clean_stack += 0.5870 * video_array[first_frame:last_frame,y : y_max,x : x_max,1]
            clean_stack += 0.1140 * video_array[first_frame:last_frame,y : y_max,x : x_max,2]
        
        ## Only cropping, no grayscaling
        else:
            clean_stack = video_array[first_frame:last_frame,y : y_max,x : x_max,:]
            
        if self.debug: print('detector.crop_and_grayscale: Final video array dimensions:',clean_stack.shape)
        return clean_stack
    
    ## Subtract background
    def subtract_background(self,video_array=None):
        '''Generate a null background image and subtract that from each frame
        ----
        Inputs:
          video_array (nd-array): clean_stack generated from crop_and_grayscale
          first_frame (int): First frame to consider for background subtraction
          last_frame (int): Last frame to consider for background subtraction
        ----
        Returns:
          spot_stack (nd-array): Background-subtracted image stack
          background (array): Array containing the pixel intensities for each x,y-coordinate'''
        if self.debug: print('detector.subtract_background')
        
        ## Setting the last frame to the end if None provided
        first_frame = self.blank_0
        last_frame = self.blank_n
                    
        ## Generating a null background image as the median pixel intensity across frames
        background = np.median(video_array[first_frame:last_frame,:,:].astype(float), axis=0).astype(int)
        if self.debug: print('detector.subtract_background: dimensions:', background.shape)
        
        ## Subtracting the null background image from each individual frame
        spot_stack = np.subtract(video_array,background)
        return spot_stack, background   


    ## Plots and views
    def view_ROI(self,image = None, border = True, x0 = 0,x1 = None, y0 = 0,y1 = None,
                 color = 'r', bin_lines = None, **kwargs):
        '''Generates image of the first frame w/a rectangle for the region of interest
        ----
        Inputs:
          image (array): Slice (single frame) of the image_stack nd-array.
          border (bool): True draws a rectangle over the region of interest
          x0 (int): Left-most coordinate
          x1 (int): Right-most coordinate
          y0 (int): Top-most coordinate (should also be the lowest y-value)
          y1 (int): Bottom-most coordinate (should also be the highest y-value)
          color (str): Color corresponding with rectangle color for region of interest
          bin_lines (bool): True if drawing lines between calculated vials
          **kwargs: Arguments for plt.imshow
        ----
        Returns:
          None -- Generates an image saved in step_3()
        '''
        if self.debug: print('detector.view_ROI')

        ## Defaults to first frame of image stack if None specified
        if self.debug: print('detector.view_ROI :: Setting frame')
        if image == None:
            image = self.image_stack[0]
        
        ## Plots the slice of nd-array
        if self.debug: print('detector.view_ROI :: Plotting image')
        plt.imshow(image,cmap=cm.Greys_r, **kwargs)
        
        ## Draws a red rectangle over the region of interest
        if self.debug: print('detector.view_ROI :: Draw ROI')
        if border:
            if x1 == None: x1 = self.image_stack[0].shape[1]
            if y1 == None: y1 = self.image_stack[0].shape[0]
            plt.hlines(y0,x0,x1, color = color, alpha = .7)
            plt.hlines(y1,x0,x1, color = color, alpha = .7)
            plt.vlines(x0,y0,y1, color = color, alpha = .7)
            plt.vlines(x1,y0,y1, color = color, alpha = .7, label='ROI')

        ## Draws box to denote where outlier trim lines are
        if self.debug: print('detector.view_ROI :: Trim outlier lines (if selected)')
        if self.trim_outliers:
            lc,rc,tc,bc = self.left_crop,self.right_crop, self.top_crop,self.bottom_crop
            plt.vlines(x0+lc,y0+bc,y0+tc,color='c',alpha=.7,linewidth=.5)
            plt.vlines(x0+rc,y0+bc,y0+tc,color='c',alpha=.7,linewidth=.5)
            plt.hlines(y0+bc,x0+lc,x0+rc,color='c',alpha=.7,linewidth=.5)
            plt.hlines(y0+tc,x0+lc,x0+rc,color='c',alpha=.7,linewidth=.5,label='Outlier trim')

        ## Draws lines between vials
        if self.debug: print('detector.view_ROI :: Drawing vial/bin lines')
        if bin_lines:
            for item in self.bin_lines[:-1]:
                plt.vlines(self.x + item, y0, y1, color = 'w', alpha = .8, linewidth = 1)
            plt.vlines(self.x + self.bin_lines[-1], y0, y1, color = 'w', alpha = .8, 
                        linewidth = 1, label='Vial boundaries')
        plt.legend()
        plt.tight_layout()
        return

    def display_images(self, cropped_converted, background, subtracted, frame=0,**kwargs):
        '''Generates a three-part figure for manipulated video frames:
          1. Cropped, converted, and grayscaled frame
          2. Null background (median pixel intensity across all indicated (blank_) frames)
          3. Result of subplots 1 - 2 (background subtracted frame)
        ----
        Inputs:
          cropped_converted (nd-array): Cropped and converted nd-array
          background (array): Null background array
          subtracted (nd-array): Background subtracted nd-array
          frame (int): Specific frame/slice of cropped_converted and subtracted
          **kwargs: Arguments for plt.imshow
        ----
        Returns:
          None -- Generates a plot saved in step_3()
          '''
        if self.debug: print('detector.displaying_images: ',end='')
        plt.figure(figsize = (6,8))
    
        ## Displaying the test frame image
        plt.subplot(311)
        if self.debug: print('| Cropped and converted |',end='')
        plt.title('Cropped and converted, frame: %s' % str(frame))
        plt.imshow(cropped_converted[frame], cmap = cm.Greys_r, **kwargs)
        plt.ylabel('Pixels')

        ## Displaying the background image
        plt.subplot(312)
        if self.debug: print(' Background image |',end='')
        plt.title('Background image')
        plt.imshow(background, cmap = cm.Greys_r, **kwargs)
        plt.ylabel('Pixels')

        ## Displaying the background subtracted, test frame image
        plt.subplot(313)
        if self.debug: print(' Subtracted background')
        plt.title('Subtracted background')
        plt.imshow(subtracted[frame], cmap = cm.Greys_r, **kwargs)
        plt.xlabel('Pixels')
        plt.ylabel('Pixels')

        plt.tight_layout()
        return

    def image_metrics(self, spots, image, metric, colorbar=False, **kwargs):
        '''Creates a plot with spot metrics placed over the video image
        ----
        Inputs:
          spots (DataFrame): DataFrame (df_big) containing the spots and their metrics
          image (array): Image background for scatter point plot
          metric (str): Spot metric to filter for
          colorbar (bool): Include a color bar legend
          **kwargs: Keyword arguments to use with plt.scatter
        ----
        Returns:
          None'''
        if self.debug: print('detector.image_metrics')

        ## Create plot
        plt.title(metric)
        plt.imshow(image, cmap = cm.Greys_r)
        plt.scatter(spots.x,spots.y, c = spots[metric], cmap = cm.coolwarm, **kwargs)
        
        ## Add in colorbar
        if colorbar: plt.colorbar()
        
        ## Format plot
        plt.ylim(self.h,0)
        plt.xlim(0,self.w)
        plt.tight_layout()
        return

    def colored_hist(self, spots, metric, bins=40, predict_threshold=False, threshold=None):
        '''Creates a colored histogram to go with image_metrics.
        ----
        Inputs:
          spots (DataFrame): DataFrame with spot metrics (df_big)
          metric (str): Spot metric to evaluate
          bins (int): Number of bins to include for histogram. Default = 40
          predict_threshold (bool): Will predict a threshold for 'signal' metric
          threshold (int): Filtering threshold
        ----
        Returns:
          None'''
        if self.debug: print('detector.colored_hist')        
        
        ## Testing threshold input value
        try: threshold = int(threshold)
        except: pass
    
        ## Assembling histogram parameters
        set_cm = plt.cm.get_cmap('coolwarm')
        n, bin_assignments, patches = plt.hist(spots[metric],bins = bins)
        bin_centers = 0.5 * (bin_assignments[:-1] + bin_assignments[1:])
        col = bin_centers - min(bin_centers)
        col /= max(col)

        ## Plotting by color
        for c, p in zip(col, patches):
            plt.setp(p, 'facecolor', set_cm(c))
    
        ## Getting height of vertical line
        y_max = np.histogram(spots[metric],bins=bins)[0].max()    

        ## Plotting vertical line for eccentricity and mass
        if metric == 'ecc': x_pos = [self.ecc_low,self.ecc_high]
        if metric == 'mass': x_pos = [self.minmass]
        if metric == 'ecc' or metric == 'mass': plt.vlines(x_pos,0,y_max,color = 'gray')
            
        ## Estimate auto-threshold
        if predict_threshold:
            _threshold = self.find_threshold(spots[metric],bins = bins)
            plt.vlines(_threshold,0,y_max,color = 'gray',label='Auto')

        ## Add in user-defined threshold vs. auto
        if isinstance(threshold, int) or isinstance(threshold, float):
            plt.vlines(threshold,0,y_max,color = 'k', label='User-defined')
            if predict_threshold:
                plt.legend(frameon=False)
        
        ## Adding y-axis labels
        plt.ylabel("Counts")
        return

    def spot_checker(self, spots, metrics=['signal'], image=None, **kwargs):
        '''Generates figure containing subplots for image_metrics and colored_hist for
          different spot metrics
        ----
        Inputs:
          spots (DataFrame): DataFrame with spot metrics (df_big)
          metrics (list): List of the metrics to include when generating the plots
          image(array): Background image for plots, default = clean_stack[0]
          **kwargs: Keyword arguments to use with plt.imshow in image_metrics function
        ----
        Returns:
          None'''
        if self.debug: print('detector.spot_checker')

        ## Setting up figure parameters
        subplot = int(str(len(metrics)) + str(2) + str(0))
        if image==None: image = self.clean_stack[0]
        count = 0
        plt.figure(figsize=(4+image.shape[1]/150,len(metrics)*2))
        
        ## Plotting each of the detector spot results over the image
        for i in range(0,len(metrics)):
            ## Defining spot metric and whether to auto-threshold
            col = metrics[i]
            if col=='signal': predict = True
            else: predict = False

            ## Drawing histogram, showing distribution of spots by metric
            count += 1
            plt.subplot(len(metrics),2,count)
            plt.title('Histogram for: %s' % col)
            plt.xlabel(col)
            self.colored_hist(spots,metric=col, bins = 40,predict_threshold=predict)
    
            ## Drawing image plot, colored by metric
            count += 1
            plt.subplot(len(metrics),2,count)
            plt.title('Spot overlay: %s' % col)
            self.image_metrics(spots,image, metric=col,**kwargs)
            plt.ylabel('Pixels')
            if col=='signal':
                plt.xlabel('Pixels')

        plt.tight_layout()
        return

    def find_spots(self, stack, diameter = 3, quiet=True,**kwargs):
        '''Locates the x,y-coordinates for all spots across frames
        ----
        Inputs:
          stack (nd-array): cropped, grayscaled, and background subtracted nd-array
          diameter (int): Estimated diameter of a spot, odds only
          quiet (bool): True silences the output
          **kwargs: Keyword arguments to use with trackpy.batch
        ----
        Returns:
          spots (DataFrame): DataFrame containing all the spots from the TrackPy output
        '''
        if self.debug: print('detector.find_spots')
        ## Check diameter
        diameter = int(diameter)
        if diameter % 2 == 0: diameter = diameter + 1
    
        ## Option to silence output
        if quiet: tp.quiet()
    
        ## Detect spots
        spots = tp.batch(stack,diameter = diameter, **kwargs)
        
        ## Sorting DataFrame
        spots = spots[spots.raw_mass > 0].sort_values(by='frame')
        return spots
                
    def particle_finder(self, invert=True, **kwargs):
        '''Finds spots and formats the resulting DataFrame. Output can be used with TrackPy.
        ----
        Inputs:
          invert (bool): True if light background, False if dark background
          **kwargs: Keyword arguments to use with trackpy.batch
        ----
        Returns:
          df (DataFrame): DataFrame containing spots and their metrics, becomes df_big'''
        if self.debug: print('detector.particle_finder')

        ## Main spot detection function
        df = self.find_spots(stack = self.spot_stack,
                             quiet=True,invert=True,
                             diameter=self.diameter,
                             minmass=self.minmass,
                             maxsize=self.maxsize)

        ## Catch if there are no spots detected
        if df.shape[0] == 0:
            print('!! Skipping video: No spots detected. Try modifying diameter, maxsize, or minmass parameters')
            raise SystemExit
        
        ## Rounding detector outputs
        df['x'] = df.x.round(2)
        df['y'] = df.y.round(2)
        df['t'] = [round(item/self.frame_rate,3) for item in df.frame]
        df['mass'] = df['mass'].astype(int)
        df['size'] = df['size'].round(3)
        df['ecc'] = df.ecc.round(3)
        df['signal'] = df.signal.round(2)
        df['raw_mass'] = df['mass'].astype(int)
        df['ep'] = df.ep.round(1)
        df['True_particle'] = np.repeat(True,df.shape[0])
        return df

    def find_threshold(self,x_array,bins=40):
        '''Auto-generates a signal threshold by finding the local minimum between two local
          maxima, or takes the average between 0 and the global maximum
        ----
        Inputs:
          x_array (list): list of all metric (signal) points to find a threshold for
          bins (int): Number of bins to search across, default = 40.
        ----
        Returns:
          threshold (int): Auto-generated threshold
        '''
        if self.debug: print('detector.find_threshold')
        
        ## Looking at the distribution of spot metrics as a histogram
        x_array = np.histogram(x_array,bins = bins)[0]

        ## Peak finding with SciPy.signal module
        peaks = find_peaks(x_array)[0]

        ## Degenerate histogram with no interior peaks -- e.g. a single tight
        ## mode from very clean, low-noise data. There is no noise/signal
        ## valley to locate, so fall back to the dominant bin instead of
        ## indexing an empty peak list (which previously raised IndexError).
        if len(peaks) == 0:
            threshold = int(np.argmax(x_array))
            if self.debug: print('                   Threshold (fallback) =',threshold)
            return threshold

        prominences = peak_prominences(x_array,peaks)
        candidates = find_peaks(x_array, prominence=np.max(prominences))[0]

        ## The prominence filter can exclude every peak (again, on near-unimodal
        ## data); fall back to the most prominent peak directly rather than
        ## indexing an empty array.
        if len(candidates) > 0:
            threshold = int(candidates[0])
        else:
            threshold = int(peaks[int(np.argmax(prominences[0]))])
        if self.debug: print('                   Threshold =',threshold)
        return threshold

    def invert_y(self,spots):
        '''Inverts spots along the y-axis. Important for converting spots indexed for an image to a plot.
        ----
        Inputs:
          spots (DataFrame): DataFrame containing a 'y' column
        ----
        Returns:
          inv_y (list): In-place list of the y-coordinates inverted'''
        if self.debug: print('detector.invert_y')
        
        ## Inverts y-axis
        inv_y = abs(spots.y - spots.y.max())
        return inv_y
    
    def get_slopes(self):
        '''Creates a dictionary with keys for vials and values for the DataFrame sliced by
        vial. It will also calculate the local linear regression for each vial and
        returns the DataFrame containing all the slopes and linear regression statistics
        for that vial.
        ----
        Inputs (Imported from the detector object):
          df (DataFrame) : DataFrame sliced by vial
          vials (int) : Number of vials in video
          window (int) : Window size to calculate local linear regression
          vial_ID (str) : Vial-specific ID, taken from the first 'vial_id_vars' of naming convention
        ----
        Returns (Exported tp the detector object):
          result (dict) : DataFrame containing the local linear regression statistics
            by vial
          vial (dict) : Dictionary of DataFrames sliced by vial, keys are vials and values
            are DataFrames'''
        if self.debug: print('detector.get_slopes')

        ## Create empty dictionaries
        self.vial,self.result = dict(),dict()

        ## Slicing DataFrame (df.filtered) by vial and assigning slices to dictionary keys (vials)
        for i in range(1,self.vials + 2):
            try:
                ## Set dict key to '1' if only 1 vial, otherwise set dict key to vial number
                if self.vials == 1 or i == self.vials + 1: self.vial[i] = self.df_filtered
                else: self.vial[i] = self.df_filtered[self.df_filtered.vial==i]
            
                ## Setting the result to the result from the local linear regression
                self.result[i] = self.local_linear_regression(self.vial[i]).iloc[0].tolist()
            
                ## Add vial_ID to the result
                if i == self.vials + 1: v = 'all'
                else: v = self._vial_label(i)

                ## Name vial_ID
                vial_ID = ['_'.join(self.vial_ID) + '_'+ str(v)]
                self.result[i] = vial_ID + self.result[i]
                
                ## Rounding results so they are more manageable and require less space.
                self.result[i][1:3] = [int(item) for item in self.result[i][1:3]]
                self.result[i][3:] = [round(item,4) for item in self.result[i][3:]]
            except:
                print('Warning:: Could not process vial %s' % i)
        return

    def get_trim_lines(self,df,edge = 'top',sensitivity=1):
        '''Calculates spacial thresholds for cropping outlier points at the edge of window
        ----
        Inputs:
          df (DataFrame): DataFrame of all points to consider
          edge (str) {'top'|'bottom','left','right'}: Which edge to trim
           sensitivity (float): How sensitive to make the thresholding
        ----
        Returns:
          crop (float): Cutoff threshold'''
        if self.debug: print('detector.get_trim_lines ::')
        for _edge in edge:
            _list,diff_list = [],[]

            # Define axis
            if edge == 'top' or edge == 'bottom': axis = 'y'
            if edge == 'left' or edge == 'right': axis = 'x'

            # Define quantile starting value
            if edge == 'left' or edge == 'bottom': quant = 0
            if edge == 'right' or edge == 'top': quant = 0.96

            ## Find quantile boundaries
            for i in range(5):
                val = df[axis].quantile(quant + i * 0.01)
                _list.append(val)

            ## Get difference between quantile boundaries
            for i in range(4):
                diff_list.append(abs(_list[i]-_list[i+1]))

            ## Calculate boundary, cutoff, and thresholds
            ## --boundary as most extreme value
            ## --cutoff as median 0-4 or 95-99 percentiles x scalar (sensitivity)
            ## --threshold as difference between boundary and cutoff
            ## --crop as value to crop points at 
        
            if 'right' in edge or 'top' in edge:
                boundary = _list[-1] # Max value
                cutoff = np.median(diff_list[:-1]) * sensitivity
                threshold = boundary - cutoff

                if cutoff < threshold: crop = boundary - cutoff
                else: crop = cutoff
                if self.debug: print('detector.get_trim_lines ::',edge,'@',crop)
                return crop

            if 'left' in edge or 'bottom' in edge:
                boundary = _list[0] # Min value
                cutoff = np.median(diff_list[1:]) * sensitivity 
                threshold = boundary + cutoff
            
                if cutoff > threshold: crop = boundary + cutoff
                else: crop = boundary
                if self.debug: print('detector.get_trim_lines ::',edge,'@',crop,'(no crop)')            
                return crop
                
    def bin_vials(self, df, vials, percentage=1,top=False, bin_lines=None):
        '''Bin spots into vials. Function takes into account all points along the x-axis, 
          and divides them into specified number of bins based on the min and max
          points in the array.
        ----
        Inputs:
          df (DataFrame):
          vials (int): Number of vials in video
        Returns:
          bin_lines (list): Binning intervals along the x-axis
          spot_assignments (pd.Series): '''
        if self.debug: print('detector.bin_vials')

        ## Bin vials, conditional for vial quantity
        if vials == 1:
            if isinstance(bin_lines, list): bin_lines = bin_lines
            else: bin_lines = [df.x.min(),df.x.max()] 
            spot_assignments = np.repeat(1,df.shape[0])
        else: ## More than 1 vial
            if isinstance(bin_lines, list): bin_lines = bin_lines
            else: bin_lines = pd.cut(df.x,vials,include_lowest=True,retbins=True)[1]

            ## Assign spots to vials
            _labels = range(1,vials+1)
            spot_assignments = pd.cut(df.x, bins=bin_lines, labels=_labels)
            spot_assignments = pd.Series(spot_assignments).astype('int')
        
            ## Checks to make sure all vials have at least one spot. Important if a middle vial is absent or vials binned incorrectly.
            counts = np.unique(spot_assignments, return_counts = True)
            for v,c in zip(counts[0], counts[1]):
                if c == 0:
                    print('Warning: vial',v,'is empty and cannot be evaluated')

        return bin_lines, spot_assignments


    def local_linear_regression(self, df, method = 'max_r'):
        '''Performs a local linear regression, using a user-defined sliding window (self.window)
        ----
        Inputs:
          df (DataFrame): DataFrame containing all formatted points (df_filtered). 'y' needs to be converted from image indexing to plot indexing
          method (str): Two options: greatest regression coefficient (max_r) or lowest error (min_err)
        ----
        Returns:
          result (DataFrame): Single-row slice of a DataFrame corresponding with the results from the local linear regression
        '''
        if self.debug: print('detector.local_linear_regression')

        ## Defining empty variables
        result_list, result = [],pd.DataFrame()
        llr_columns = ['first_frame','last_frame','slope','intercept','r','pval','err']#,'count_llr','count_all']
        
        _count_all = np.median(df.groupby('frame').frame.count())
        
        ## Iterating through the window
        frames = (self.crop_n - self.crop_0) - self.window
        for i in range(frames):
            ## Defining search parameters for each iteration
            start, stop = int(i),int(i+self.window)
            df_window  = df[(df.frame >= start) & (df.frame <= stop)]

            ## Testing if there are enough frames in the slice
            if df_window.groupby('frame').y.count().min() == 0:
                print('Issue with number of frames with flies detected vs. window size')
                print(i, i + self.window, len(df_window.frame.unique()))
                continue
                        
            ## Performing local linear regression
            try: 
                ## Grouping points by frame
                _frame = df_window.groupby('frame').frame.mean()
                _pos  = df_window.groupby('frame').y.mean()
                _count_llr = np.median(df_window.groupby('frame').frame.count())

                ## Performing linear regression on subset and formatting output to list
                _result = linregress(_frame,_pos)
                _result = [start,stop] + np.hstack(_result).tolist() #+ [_count_llr,_count_all]

                ## If slope is not significantly different from 0, then set slope = 0
                if _result[-2] >= 0.05: _result[2] = 0

            ## Have row of NaN if unable to process
            except: _result = [start,stop] + [np.nan,np.nan,np.nan,np.nan]
#             except: _result = [start,stop] + [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

            ## Add results list to a list of lists
            result_list.append(_result)
        
        ## Assembles the list of lists into a DataFrame
        result = pd.DataFrame(data=result_list,columns=llr_columns)

        ## Filtering method
        if method == 'max_r': result = result[result.r == result.r.max()]
        elif method == 'min_err': result = result[result.err == result.err.min()]
        else: print("Unrecognized method, chose either 'max_r' to select the window with the greatest R or 'min_err' to select the window with the lowest error")
        return result

    # ---------- FNG detection helpers ----------
    def _height_traces(self):
        """
        Build per-vial mean y-position vs frame from self.df_filtered.
        Returns a DataFrame indexed by frame with one column per vial.
        """
        df = getattr(self, 'df_filtered', None)
        if df is None or df.empty:
            return pd.DataFrame()
        H = (df.groupby(['frame','vial']).y.mean()
               .unstack('vial')
               .sort_index())
        return H  # rows=frame, cols=vial (1..self.vials)

    def _detect_fng_series(self, series,
                           smooth_window=None,
                           climb_thresh=None,
                           fall_thresh=None,
                           min_gap=None):
        """
        Detect 'climb–then–fall' events in one vial's height trace.
        All thresholds operate on 0–1 normalized signal.
        Returns a list of dicts with the peak frame and rise/drop magnitudes.
        """
        # Defaults (safe if not present in cfg)
        smooth_window = smooth_window if smooth_window is not None else getattr(self, 'fng_smooth_window', 5)
        climb_thresh  = climb_thresh  if climb_thresh  is not None else getattr(self, 'fng_climb_thresh', 0.10)
        fall_thresh   = fall_thresh   if fall_thresh   is not None else getattr(self, 'fng_fall_thresh', 0.10)
        min_gap       = min_gap       if min_gap       is not None else getattr(self, 'fng_min_gap', 5)

        # Smooth and keep length (centered rolling)
        s = pd.Series(series).rolling(window=max(1, int(smooth_window)), center=True).mean().bfill().ffill()

        # Normalize 0..1 so thresholds are comparable across rigs
        rng = (s.max() - s.min())
        if not pd.notnull(rng) or rng == 0:
            return []
        n = (s - s.min()) / rng
        nv = n.values

        # Candidate peaks (tops of climbs)
        peaks, _ = find_peaks(nv)
        sv = s.values  # raw (smoothed) pixel values, same positional alignment as nv

        events = []
        last_event_frame = -10**9
        w = max(1, int(smooth_window)) * 2

        for p in peaks:
            if p - last_event_frame < min_gap:
                continue

            # Rise: measure from end of previous fall (or series start) to current peak.
            # This captures the full climb regardless of how gradually the fly ascended.
            left_start = max(0, last_event_frame + 1)
            left = nv[left_start:p+1]
            right_slice = nv[p:min(len(nv), p + w)]

            if len(left) < 2 or len(right_slice) < 2:
                continue

            left_min  = left.min()
            right_min = right_slice.min()

            rise = nv[p] - left_min
            drop = nv[p] - right_min

            if rise >= climb_thresh and drop >= fall_thresh:
                # Map positional indices back to original frame numbers
                n_valid = getattr(self, 'n_frames', None)
                _clamp = (lambda f: max(0, min(int(f), n_valid - 1))
                          if n_valid and n_valid > 0 else int(f))

                frame_peak = _clamp(n.index[p])
                p_end = p + int(np.argmin(right_slice))
                frame_fall_end = _clamp(n.index[p_end])
                drop_px = float(sv[p] - sv[p_end])
                drop_cm = drop_px / float(getattr(self, 'pixel_to_cm', 1.0))

                # Find fall start: first frame after peak where signal drops below peak value
                search_end = min(len(nv), p + w)
                if p + 1 < search_end:
                    below_idx = np.where(nv[p+1:search_end] < nv[p])[0]
                    p_fall_start = p + 1 + int(below_idx[0]) if len(below_idx) > 0 else p + 1
                else:
                    p_fall_start = p
                frame_fall_start = _clamp(n.index[min(p_fall_start, len(n) - 1)])

                # Fall duration
                frame_rate = float(getattr(self, 'frame_rate', 1.0))
                fall_duration_frames = frame_fall_end - frame_fall_start
                fall_duration_sec = round(fall_duration_frames / frame_rate, 4)

                # Recovery detection: first run of 3+ consecutive frame-over-frame
                # increases in the smoothed signal after the fall end.
                # Search is bounded by the next detected peak (or end of series).
                next_peaks = peaks[peaks > p]
                recovery_bound = int(next_peaks[0]) if len(next_peaks) > 0 else len(nv)

                consecutive = 0
                run_start_pos = None
                frame_recovery_start = None
                for i in range(p_end + 1, recovery_bound):
                    if sv[i] > sv[i - 1]:
                        if consecutive == 0:
                            run_start_pos = i - 1
                        consecutive += 1
                        if consecutive >= 3:
                            frame_recovery_start = int(n.index[run_start_pos])
                            break
                    else:
                        consecutive = 0
                        run_start_pos = None

                if frame_recovery_start is not None:
                    recovery_duration_sec = round(
                        (frame_recovery_start - frame_fall_end) / frame_rate, 4)
                else:
                    recovery_duration_sec = float('nan')

                events.append({
                    'frame_peak': frame_peak,
                    'frame_fall_start': frame_fall_start,
                    'frame_fall_end': frame_fall_end,
                    'fall_duration_frames': fall_duration_frames,
                    'fall_duration_sec': fall_duration_sec,
                    'rise_norm': float(rise),
                    'drop_norm': float(drop),
                    'drop_px': round(drop_px, 2),
                    'drop_cm': round(drop_cm, 4),
                    'recovery_duration_sec': recovery_duration_sec,
                })
                last_event_frame = p_end

        return events

    def compute_fng(self):
        """
        Compute FNG events per vial from df_filtered and save:
          - self.df_fng: event-level table
          - self.df_fng_counts: counts per vial
        Writes <video>.fng.csv alongside the other outputs.
        """
        if not getattr(self, 'fng_enabled', True):
            self.df_fng = pd.DataFrame()
            self.df_fng_counts = pd.DataFrame({'vial': range(1, getattr(self, 'vials', 0)+1),
                                               'fng_count': 0})
            return

        FNG_COLUMNS = ['vial','event_idx','frame_peak','frame_fall_start','frame_fall_end',
                       'fall_duration_frames','fall_duration_sec',
                       'rise_norm','drop_norm','drop_px','drop_cm','recovery_duration_sec']

        H = self._height_traces()
        if H.empty:
            self.df_fng = pd.DataFrame(columns=FNG_COLUMNS)
            self.df_fng_counts = pd.DataFrame({'vial': range(1, getattr(self, 'vials', 0)+1),
                                               'fng_count': 0})
            # still emit a headered csv for consistency
            path_fng = self.name_nosuffix + '.fng.csv'
            self.df_fng.to_csv(path_fng, index=False)
            print('                --> Saved:', path_fng.split('/')[-1])
            return

        records = []
        for vial in H.columns:
            series = H[vial]
            evs = self._detect_fng_series(
                series,
                smooth_window=getattr(self, 'fng_smooth_window', 5),
                climb_thresh=getattr(self, 'fng_climb_thresh', 0.10),
                fall_thresh=getattr(self, 'fng_fall_thresh', 0.10),
                min_gap=getattr(self, 'fng_min_gap', 5),
            )
            for idx, ev in enumerate(evs, start=1):
                records.append({
                    'vial': int(vial),
                    'event_idx': idx,
                    'frame_peak': ev['frame_peak'],
                    'frame_fall_start': ev['frame_fall_start'],
                    'frame_fall_end': ev['frame_fall_end'],
                    'fall_duration_frames': ev['fall_duration_frames'],
                    'fall_duration_sec': ev['fall_duration_sec'],
                    'rise_norm': round(ev['rise_norm'], 4),
                    'drop_norm': round(ev['drop_norm'], 4),
                    'drop_px': ev['drop_px'],
                    'drop_cm': ev['drop_cm'],
                    'recovery_duration_sec': ev['recovery_duration_sec'],
                })

        self.df_fng = (pd.DataFrame.from_records(records, columns=FNG_COLUMNS)
                       .sort_values('frame_peak')
                       .reset_index(drop=True))

        counts = (self.df_fng.groupby('vial').size()
                    .rename('fng_count')
                    .reindex(range(1, self.vials+1), fill_value=0)
                    .reset_index())
        self.df_fng_counts = counts

        # Save per-event csv (vial column relabeled to physical IDs if configured)
        path_fng = self.name_nosuffix + '.fng.csv'
        if self.df_fng.empty:
            pd.DataFrame(columns=FNG_COLUMNS).to_csv(path_fng, index=False)
        else:
            self._relabel_vial_col(self.df_fng).to_csv(path_fng, index=False)
        print('                --> Saved:', path_fng.split('/')[-1])

    def compute_tortuosity(self):
        """
        Compute per-fly, per-bout tortuosity metrics and save the two
        tortuosity CSV files. Runs only in individual mode (gated by the
        step_5 hook).

        Bout segmentation is INDEPENDENT of FNG events (Fix #2): this method
        does NOT call _detect_fng_series, does NOT read self.df_fng /
        self.fng_events, and does NOT use any FNG-derived data. Bouts are
        detected per particle from vertical climbing velocity.

        Configuration (all read via getattr with safe defaults so existing
        .cfg files keep working unchanged):
          tortuosity_velocity_threshold    mm/s threshold for a climbing step
          tortuosity_bout_min_frames       minimum bout length, in frames
          tortuosity_bout_min_displacement minimum bout net vertical
                                           displacement, in mm

        Writes:
          <video>.tortuosity_bouts.csv    one row per (vial, particle,
                                          bout_idx): tortuosity, straightness,
                                          vertical_efficiency, mean turning
                                          angle (rad), plus bout duration and
                                          path length in mm.
          <video>.tortuosity_particle.csv one row per (vial, particle):
                                          n_bouts and
                                          median_vertical_efficiency.

        See scripts/tortuosity.py for metric definitions and the
        velocity-threshold bout-segmentation algorithm.
        """
        if self.debug: print('detector.compute_tortuosity')

        df = getattr(self, 'df_filtered', None)

        # Config (Fix #3 allowlist keys, defensive getattr with safe defaults
        # so any .cfg file that omits the tortuosity_* keys keeps working).
        smoothing_window      = getattr(self, 'tortuosity_smoothing_window', 5)
        velocity_threshold    = getattr(self, 'tortuosity_velocity_threshold', 1.0)
        bout_min_frames       = getattr(self, 'tortuosity_bout_min_frames', 10)
        bout_min_displacement = getattr(self, 'tortuosity_bout_min_displacement', 2.0)
        pixel_to_cm           = getattr(self, 'pixel_to_cm', 1.0)
        frame_rate            = getattr(self, 'frame_rate', 1.0)

        if df is None or df.empty or 'particle' not in df.columns:
            print('   No linked tracks; skipping tortuosity computation')
            self.df_tortuosity_bouts = pd.DataFrame(
                columns=_tortuosity.TORTUOSITY_BOUT_COLUMNS)
        else:
            print('-- [ Tortuosity ] Computing per-fly bout metrics '
                  '(velocity threshold = %.3f mm/s)' % velocity_threshold)
            self.df_tortuosity_bouts = _tortuosity.compute_tortuosity_table(
                df,
                pixel_to_cm=pixel_to_cm,
                frame_rate=frame_rate,
                velocity_threshold=velocity_threshold,
                bout_min_frames=bout_min_frames,
                bout_min_displacement=bout_min_displacement,
                smoothing_window=smoothing_window,
            )
            print('   %d bout-particle row(s) computed'
                  % len(self.df_tortuosity_bouts))

        self.df_tortuosity_particle = _tortuosity.compute_particle_table(
            self.df_tortuosity_bouts)

        # Backward-compat alias used by earlier callers / tests.
        self.df_tortuosity = self.df_tortuosity_bouts

        path_bouts = self.name_nosuffix + '.tortuosity_bouts.csv'
        path_particle = self.name_nosuffix + '.tortuosity_particle.csv'
        self._relabel_vial_col(self.df_tortuosity_bouts).to_csv(path_bouts, index=False)
        self._relabel_vial_col(self.df_tortuosity_particle).to_csv(path_particle, index=False)
        print('                --> Saved:', path_bouts.split('/')[-1])
        print('                --> Saved:', path_particle.split('/')[-1])
        return

    def link_trajectories(self):
        """
        Link per-frame detections into per-fly trajectories, one vial at a time.

        Runs only when analysis_mode == 'individual' (the step_5 hook gates this
        call). Operates on self.df_filtered -- whose frame/vial/x/y columns are
        already populated by step_4/step_5 -- assigns a globally unique
        'particle' ID to every surviving detection, drops fragmentary tracks via
        tp.filter_stubs, and writes <video>.tracks.csv alongside the other
        outputs.

        Modeled structurally on compute_fng(): config is pulled with defensive
        getattr() so existing .cfg files lacking the link_* keys keep working,
        each vial is processed independently, results are stored back on the
        detector, and a CSV is saved. Cohort-mode analysis never reaches here,
        so existing FNG detection output is unaffected.
        ----
        Inputs:
          None -- reads self.df_filtered and the link_* configuration keys
        ----
        Returns:
          None -- self.df_filtered gains a 'particle' column; <video>.tracks.csv
                  is written
        """
        if self.debug: print('detector.link_trajectories')

        ## Defensive config access -- mirrors the fng_* default pattern
        search_range = getattr(self, 'link_search_range', 15)
        memory       = getattr(self, 'link_memory', 3)
        predictor    = getattr(self, 'link_predictor', 'nearest_velocity')
        min_len      = getattr(self, 'link_min_track_length', 5)

        print('-- [ Linking ] Linking per-fly trajectories (predictor: %s)' % predictor)

        df = getattr(self, 'df_filtered', None)
        if df is None or df.empty:
            print('   No detections to link; skipping trajectory linking')
            return

        TRACK_COLUMNS = ['particle', 'frame', 't', 'vial', 'x', 'y']
        naming_cols = list(getattr(self, 'file_details', {}).keys())

        linked_frames = []
        particle_offset = 0
        total_stubs = 0
        vials_processed = 0

        ## Link each vial INDEPENDENTLY so particle IDs never swap across vials.
        ## The 'vial' column was assigned in step_4.
        for vial, group in df.groupby('vial'):
            g = group.copy()

            ## Branch on the configured predictor
            if predictor == 'nearest_velocity':
                pred = tp.predict.NearestVelocityPredict()
                linked = pred.link_df(g, search_range=search_range, memory=memory)
            elif predictor == 'none':
                linked = tp.link_df(g, search_range=search_range, memory=memory)
            else:
                # TODO: unrecognized link_predictor value -- spec allows only
                # 'nearest_velocity' or 'none'; falling back to plain linking.
                print('   !! Unknown link_predictor (%s); using plain tp.link_df' % predictor)
                linked = tp.link_df(g, search_range=search_range, memory=memory)

            ## Drop fragmentary tracks shorter than link_min_track_length frames
            n_before = linked['particle'].nunique()
            linked = tp.filter_stubs(linked, threshold=min_len).reset_index(drop=True)
            n_after = linked['particle'].nunique()
            n_stubs = n_before - n_after
            total_stubs += n_stubs
            vials_processed += 1

            if linked.empty:
                print('   vial %s: 0 tracks survived (%s stub track(s) filtered)'
                      % (vial, n_stubs))
                continue

            ## Offset particle IDs so each vial occupies a disjoint ID range
            linked['particle'] = linked['particle'].astype(int) + particle_offset
            particle_offset = int(linked['particle'].max()) + 1

            n_particles = linked['particle'].nunique()
            mean_len = linked.groupby('particle').size().mean()
            print('   vial %s: %s unique particle(s), %s stub track(s) filtered, '
                  'mean track length %.1f frames'
                  % (vial, n_particles, n_stubs, mean_len))

            linked_frames.append(linked)

        ## Concatenate per-vial results back into df_filtered (now with 'particle')
        if linked_frames:
            self.df_filtered = pd.concat(linked_frames, ignore_index=True)
        else:
            self.df_filtered = df.iloc[0:0].copy()
            self.df_filtered['particle'] = pd.Series(dtype=int)

        print('-- [ Linking ] %s vial(s) processed, %s total stub track(s) filtered'
              % (vials_processed, total_stubs))

        ## Write *.tracks.csv mirroring the *.filtered.csv naming convention
        track_cols = [c for c in TRACK_COLUMNS + naming_cols
                      if c in self.df_filtered.columns]
        path_tracks = self.name_nosuffix + '.tracks.csv'
        self._relabel_vial_col(self.df_filtered[track_cols]).to_csv(path_tracks, index=False)
        print('                --> Saved:', path_tracks.split('/')[-1])
        return

    def step_1(self, gui = False, grayscale = True):
        '''Crops and formats the video, previously loaded during detector initialization.
        ----
        Inputs:
          gui (bool): True creates plots for detector optimization
          grayscale (bool): True converts video array to grayscale
        ----
        Returns:
          None'''
        print('-- [ Step 1  ] Cleaning and format image stack')

        # Release large arrays from any prior run before allocating new ones so
        # the old and new stacks do not coexist during the crop_and_grayscale call.
        for _attr in ('clean_stack', 'spot_stack', 'background'):
            if hasattr(self, _attr):
                delattr(self, _attr)
        gc.collect()

        x,y = self.x,self.y
        x_max, y_max = int(x + self.w),int(y + self.h)
        stack = self.image_stack

        ## Confirm frame ranges
        self.check_variable_formats()
        if self.blank_0 < self.crop_0:
            self.blank_0 = self.crop_0
        if self.blank_n > self.crop_n:
            self.blank_n = self.crop_n
        
        if grayscale:
            if self.debug: print('detector.step_1 cropped and grayscale: grayscale image')
            self.clean_stack = self.crop_and_grayscale(stack,
                         y=y, y_max=y_max,
                         x=x, x_max=x_max,
                         first_frame=self.crop_0, last_frame=self.crop_n)
        else:
            if self.debug: print('detector.step_1 cropped and grayscale: no color image')
            self.clean_stack = self.crop_and_grayscale(stack,
                         y=y, y_max=y_max,
                         x=x, x_max=x_max,
                         first_frame=self.crop_0, last_frame=self.crop_n, grayscale=False)                        

        if self.debug: print('detector.step_1 cropped and grayscale dimensions: ', self.clean_stack.shape)

        ## Subtracts background to generate null background image and spot stack
        self.spot_stack,self.background = self.subtract_background(video_array=self.clean_stack)
        if self.debug: print('detector.step_1 spot_stack and null background created')
        return


    def step_2(self):
        '''Performs spot detection and manipulates the resulting DataFrames'''
        print('-- [ Step 2  ] Identifying spots')

        ## Particle detection step
        self.df_big = self.particle_finder(minmass=self.minmass,diameter=self.diameter,
                                            maxsize=self.maxsize, invert=True)
        if self.debug: print('                   Identified %s spots' % self.df_big.shape[0])
        return


    def step_3(self, gui = False):
        '''Visualizes spot metrics
        ----
        Inputs:
          gui (bool): True creates plots for detector optimization
        ----
        Returns:
          None
        '''
        print('-- [ Step 3  ] Visualize spot metrics ::',gui)
        if gui:        
            ## Visualizes spot metrics on plot with accompanying color-coded histogram
            self.spot_checker(self.df_big,metrics=['ecc','mass','signal'], alpha=.1)
            
            plot_spot_check = self.name_nosuffix + '.spot_check.png'
            plt.savefig(plot_spot_check,dpi=200)
            print('                --> Saved:',plot_spot_check.split('/')[-1])
            plt.close()
        
            ## Creating image plot with rectangle superimposed over first frame
            plt.figure()
            self.display_images(self.clean_stack,self.background,self.spot_stack,frame=20)
            plt.tight_layout()
            plot_name = self.name_nosuffix + '.processed.png'
            plt.savefig(plot_name, dpi=100)
            plt.close()
            print('                --> Saved:',plot_name.split('/')[-1])
        return

    def step_4(self):
        '''Filters and processes data detected points'''

        ## Assigning spots a True/False status based on ecc/eccentricity (circularity)
        def ecc_filter(x,low=0,high=1):
            '''Simple function for making a vector True/False depending on an upper and lower bound
            ----
            Inputs:
              x (numeric): value to perform operation on
              low (numeric): Lower bound
              high (numeric): Upper bound
            Returns:
              (bool): True/False'''
            if x >= low and x <= high: return True
            else: return False

        print('-- [ Step 4a ]   - Setting spot threshold')        
        ## Auto-detecting threshold
        if self.threshold == 'auto': self.threshold = self.find_threshold(self.df_big.signal)

        print('-- [ Step 4b ]   - Filtering by signal threshold') 
        ## Assigning spots a True/False status based on signal threshold
        self.df_big['True_particle'] = [x >= self.threshold for x in self.df_big.signal]

        t_or_f = np.unique(self.df_big.True_particle, return_counts=True)
        if self.debug: print('                   True (%s) and False (%s) spots' % (t_or_f[0],t_or_f[1]))
        
        print('-- [ Step 4c ]   - Filtering by eccentricity/circularity') 
        self.df_big.loc[self.df_big.True_particle == True,'True_particle'] = self.df_big[self.df_big.True_particle==True].ecc.map(lambda x: ecc_filter(x,low=self.ecc_low,high=self.ecc_high))
        t_or_f = np.unique(self.df_big.True_particle, return_counts=True)
        if self.debug: print('                   True (%s) and False (%s) spots'%(t_or_f[0],t_or_f[1]))
        
        ## Checking to confirm DataFrame is not empty after filtering
        if self.df_big[self.df_big.True_particle].shape[0] == 0:
            print('\n\n!! No spots post-filtering, check detector and background subtraction settings for proper optimization')
            raise SystemExit

        ## Pruning errant points on periphery if outliers
        print('-- [ Step 4d ]   - Trimming outliers (if indicated)')
        if self.trim_outliers:
            self.left_crop = self.get_trim_lines(self.df_big,edge='left',sensitivity = self.outlier_LR)
            self.right_crop = self.get_trim_lines(self.df_big,edge='right',sensitivity = self.outlier_LR)
            self.top_crop = self.get_trim_lines(self.df_big,edge='top',sensitivity = self.outlier_TB)
            self.bottom_crop = self.get_trim_lines(self.df_big,edge='bottom',sensitivity = self.outlier_TB)
            
            self.df_big = self.df_big[(self.df_big.x >= self.left_crop) & (self.df_big.x <= self.right_crop) &
                                        (self.df_big.y <= self.top_crop) & (self.df_big.y >= self.bottom_crop)]
        
        ## Assigning spots to vials, 0 if False AND outside of the True point range
        print('-- [ Step 4e ]   - Assigning spots to vials')
        self.bin_lines, self.df_big.loc[self.df_big['True_particle'],'vial'] = self.bin_vials(self.df_big[self.df_big.True_particle],vials = self.vials)
        
        ########################################
        ## Beginning of publication insert
        if publication:
            df=self.df_big
            self.bin_lines = self.bin_vials(df[df.y > 120], vials= self.vials)[0]
            df['vial'] = np.repeat(0,df.shape[0])
            vial_assignments = self.bin_vials(df, vials = self.vials, bin_lines = self.bin_lines)[1]
            df.loc[(df.x >= self.bin_lines[0]) & (df.x <= self.bin_lines[-1]),'vial'] = vial_assignments
            self.df_big = df
        ## End of publication insert
        ########################################
        
        self.df_big.loc[self.df_big['True_particle']==False,'vial'] = 0

        print('-- [ Step 4f ]   - Saving raw data file')
        ## Saving the TrackPy results, plus filter and vial notations        
        self.df_big.to_csv(self.path_data, index=None)
        print('                --> Saved:',self.path_data.split('/')[-1])

        return

    def step_5(self):
        '''Calculates local linear regressions'''
        print('-- [ step 5  ] Setting up DataFrames for local linear regression')

        ## Filtering spots and pruning unnecessary columns
        self.df_filtered = self.df_big[(self.df_big.True_particle) & (self.df_big.vial != 0)]
        if self.debug: print('self.df_filtered.shape:',self.df_filtered.shape)
        self.df_filtered = self.df_filtered.drop(['ecc','signal','ep','raw_mass','mass','size','True_particle'],axis=1)
        self.df_filtered = self.df_filtered.sort_values(by=['vial','frame','y','x'])    

        ## Adding experimental details to DataFrame
        self.specify_paths_details(self.video_file)
       
       ## Filling in experimental details to DataFrame
        for item in self.file_details.keys():
            self.df_filtered[item] = np.repeat(self.file_details[item],self.df_filtered.shape[0])        
        
        ## Invert y-axis -- images indexed upper left to lower right but converting because plots got left left to upper right
        self.df_filtered['y'] = self.invert_y(self.df_filtered)
        self.df_filtered['y'] = self.df_filtered.y.round(2)
        
        ## Convert vial assignments from float to int
        self.df_filtered['vial'] = self.df_filtered['vial'].astype('int')

        #---- Per-fly trajectory linking (individual mode only) ----
        if getattr(self, 'analysis_mode', 'cohort') == 'individual':
            self.link_trajectories()

        #----FNG detection (per vial) ----
        self.compute_fng()

        #---- Per-fly tortuosity metrics (individual mode + tortuosity_enabled) ----
        if (getattr(self, 'analysis_mode', 'cohort') == 'individual'
                and getattr(self, 'tortuosity_enabled', True)):
            self.compute_tortuosity()

        ## Save the filtered DataFrame
        path_filtered = self.name_nosuffix+'.filtered.csv'
        self.df_filtered.to_csv(self.path_filtered, index=False)
        print('                --> Saved:',self.path_filtered.split('/')[-1])
        return
        
    def step_6(self,gui=False):
        '''Creating diagnostic plots to visualize spots at beginning & end of most linear 
          section, throughout the video, and a vertical velocity plot for each vial.
        ----
        Inputs:
          gui (bool): True creates plots for detector optimization
        ----
        Returns:
          None'''
        print('-- [ step 6a ] Visualize spot metrics ::',gui)
        
        ## Check window size is not greater than video length
        video_length = self.crop_n - self.crop_0
        if self.window > video_length:
            print('!! Issue with window size > video length: was %s, now %s' % (self.window, video_length-1))
            self.window = video_length - 1
            
        if gui:        

            ## Visualizes the region of interest
            plt.figure()
            self.view_ROI(border = True,
                            x0 = self.x, x1 = self.x + self.w,
                            y0 = self.y, y1 = self.y + self.h,
                            bin_lines = True)
        
            plot_roi = self.name_nosuffix + '.ROI.png'
            plt.savefig(plot_roi,dpi=100)
            print('                --> Saved:',plot_roi.split('/')[-1])
            plt.close()
            
        print('-- [ step 6b ] Creating diagnostic plot file')
        ## Set up plots
        plt.figure(figsize=(10,8))
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

        ## Finding the frames that flank the most linear portion 
        ##    of the y vs. t curve for all points, not just by vials
        if self.debug: print('-- [ step 6b ] Plotting data: Re-running local linear regression on all')
        _result = self.local_linear_regression(self.df_filtered)
        begin = _result.iloc[0].first_frame.astype(int)
        end = _result.iloc[0].last_frame.astype(int)
        
        ## For future release
#         min_R = _result.iloc[0].r_value ##

        ## Need only True spots, but not inverted -- df_big
        spots = self.df_big[self.df_big.True_particle]

        ## Creating the diagnostic plot
        print('-- [ step 6b1] Plotting image plots with overlaying points')        
        if self.debug: print("-- [ step 6b1] Plotting data: Plot 1 - Frame %s" % begin)
        self.image_plot(df = spots,frame = begin, ax=ax1)

        if self.debug: print("-- [ step 6b1] Plotting data: Plot 2 - Frame %s" % end)
        self.image_plot(df = spots, ax=ax2, frame = int(str(end)))

        if self.debug: print("-- [ step 6b1] Plotting data: Plot 3 - Frame ALL")
        self.image_plot(df = spots,ax=ax4, frame = None)

        print('-- [ step 6b2] Performing local linear regression')
        self.get_slopes()
        
        print('-- [ step 6b3] Plotting local linear regression results')
        self.loclin_plot(ax=ax3)

        ## Saving diagnostic plot
        fig.tight_layout()
        plt.savefig(self.path_diagnostic,dpi=300, transparent = True)
#         plt.savefig(self.path_project + '/diagnostic_plots/' + self.name + '.diagnostic.png',dpi=100, transparent = True)        

        ## Future release
#         if min_R < self.review_R:
#             plt.savefig(self.path_review_diagnostic,dpi=100, transparent = True)
            
        plt.close()
        plt.close()
        print('                --> Saved:',self.path_diagnostic.split('/')[-1])
        return
        
    def step_7(self):
        '''Writing the video's slope file'''
        print('-- [ step 7  ] Setting up slopes file')
        slope_columns = ['vial_ID','first_frame','last_frame','slope','intercept','r_value','p_value','std_err']#,'count_llr','count_all']
        
        ## Converting dictionary of local linear regressions into a DataFrame
        self.df_slopes = pd.DataFrame.from_dict(self.result,orient='index',columns = slope_columns)

        ## Applying conversion factor if indicated, if not it will just be '1'
        self.df_slopes['slope'] = self.df_slopes.slope.transform(lambda x: x * (self.conversion_factor)).round(4)
        
        ## Adding in experimental details from naming convention into the slopes DataFrame
        for item in self.file_details.keys():
            self.df_slopes[item] = np.repeat(self.file_details[item],self.df_slopes.shape[0])

        ## Specifying column names
        slope_columns = [item for item in self.file_details.keys()] + slope_columns
        self.df_slopes = self.df_slopes[slope_columns]
    
        ## Saving slope file
        self.df_slopes.to_csv(self.path_slope,index=False)
        print('                --> Saved: %s \n' % self.path_slope.split('/')[-1])
        plt.close('all')
        
        print(self.df_slopes[['vial_ID','slope','r_value']])
        print('\n')
        return
        
    def image_plot(self,df,frame=None,ax=None,ylim=[0,1000]):
        '''Image subplot for the diagnostic plot
        ----
        Inputs:
          df (DataFrame): DataFrame containing all spots
          frame (int): Frame to slice df
          ax (int): plot coordinate
          ylim (2-item list): y-limits
        ----
        Returns:
          ax (object): matplotlib object containing plot'''
        if self.debug: print('detector.image_plot')
        
        ## Get frame number
        try: frame = int(frame)
        except: frame = None

        ## Assign plotting parameters depending on which frame(s)
        if isinstance(frame, (int, np.integer)) and frame in df.frame.unique():
            df = df[(df.frame == frame)]
            alpha = .25
            title  = "Frame: %s" % frame
            ax.set_title(title)
        elif frame == None:
            frame = 0
            alpha = 0.01
            title = 'All x,y-points throughout video'
            ax.set_title(title)
        elif isinstance(frame, (int, np.integer)) and frame not in df.frame.unique():
            print('Error: input frame is not accounted for in current DataFrame')
        else:
            print("Chose a frame value in integer form, or 'None'")

        ## Issue with plotting if last frame in stack
        if frame == self.n_frames: image = self.clean_stack[frame-1]
        else: image = self.clean_stack[frame]

        ## Plotting image
        ax.imshow(image,cmap=cm.Greys_r,origin='upper')
        ax.set_ylim(self.h,0)
        ax.set_xlim(0,self.w)
        
        ## Plotting vertical bin lines
        ax.vlines(self.bin_lines,0,image.shape[0],alpha = .3)  
        
        ## Coloring spots by vial
        df = df.sort_values(by='vial')
        df = df[df.vial != 0]
        if self.vials >= 1:
            ax.scatter(df.x, df.y, 
                        s = 30, 
                        alpha = alpha,
                        c = df['vial'], 
                        cmap=self.vial_color_map)

        ## Getting rid of axis labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        return ax
    
    def loclin_plot(self,ax = None):
        '''Local linear regression plot: mean y-position vs. frame or time. Colored by vial
          and bolded for the most linear section
        ----
        Inputs:
          ax (int): plot coordinate
        ----
        Returns:
          None'''
        if self.debug: print('detector.loclin_plot')
        
        def two_plot(df,vial,label,first,last,ax=None):
            '''Adds the bolded flair to the local linear regression plot'''
            if self.debug: print('detector.two_plot')
            
            ## All points
            x = df.groupby('frame').frame.mean()
            y = df.groupby('frame').y.mean()
            
            ## Convert to cm / sec
            if self.convert_to_cm_sec:
                x = x / self.frame_rate
                y = y / self.pixel_to_cm
            ax.plot(x,y, 
                    alpha = .35, 
                    color = self.color_list[vial-1],
                    label='') 

            ## Only points in the most linear segment
            df = df[(df.frame >= first) & (df.frame <= last)]
            x = df.groupby('frame').frame.mean()
            y = df.groupby('frame').y.mean()
            
            ## Convert to cm / sec
            if self.convert_to_cm_sec:
                x = x / self.frame_rate
                y = y / self.pixel_to_cm
            ax.plot(x,y,
                    color = self.color_list[vial-1],
                    label = label)
            return
    
        ## Plotting multiple vials' data
        if len(self.result) > 1:
            for V in range(1,self.vials+1):
                l = 'Vial '+str(V)
                c = self.color_list[V-1]
                _details = self.result[V]
                frames = _details[1],_details[2]
                two_plot(self.vial[V],
                         vial = V,
                         label = l,
                         first = _details[1],
                         last  = _details[2],
                         ax=ax)

        ## Add on labels and legends
        ax.legend(loc=2, frameon=False, fontsize='x-small')
        label_y,label_x = 'pixels','Frames'
        if self.convert_to_cm_sec: 
            label_x = 'Seconds'
            label_y = 'cm'
        title = 'Cohort climbing kinematics'
        label_y = 'Mean y-position (%s)' % label_y
        ax.set_ylim(0,ax.get_ylim()[1])
        ax.set(title=title,xlabel = label_x,ylabel=label_y)
        
        if ax == None: return
        else: return ax

    ## Parameter testing is only used in the GUI
    def parameter_testing(self, variables, axes):
        '''Parameter testing in the GUI and done separately to account for plots with wx'''
        if self.debug: print('detector.parameter_testing')

        ## Close any pyplot figures left open by a prior run before allocating new ones.
        plt.close('all')

        ## Running through the first few steps
        self.load_for_gui(variables)

        ## Load in video
        self.step_1(gui=True) # Crop and convert video

        ## Detect spots
        self.step_2()
        
        ## First optimization plot
        self.step_3(gui=True)
        
        ## Filters DataFrame of detected spots
        self.step_4() 

        ## Executing final steps
        self.step_5()
        self.step_6(gui=True)
        self.step_7()
        
        #### Working through the GUI plots        
        ## Setting plot (upper left) for background image
        if self.debug: print('detector.parameter_testing: Subplot 0: Background image')
        axes[0].set_title("Background Image")
        axes[0].imshow(self.background,cmap=cm.Greys_r)
        axes[0].set_xlim(0,self.w)
        axes[0].set_ylim(self.h,0)
        axes[0].scatter([0,self.w],[0,self.h],alpha=0,marker='.')
        
        ## Slice df_big into true vs. false spots
        if self.debug: print('detector.parameter_testing: Slicing DataFrames')
        spots_false = self.df_big[~self.df_big['True_particle']]
        spots_true = self.df_big[self.df_big['True_particle']]
        
        ## Binning and coloring spots
#         bin_lines,spots_true['vial'] = self.bin_vials(spots_true,vials = self.vials)
#         spots_true = spots_true[(spots_true.x >= self.bin_lines.min()) & (spots_true.x <= self.bin_lines.max())]

        spots_true['vial'] = np.repeat(0,spots_true.shape[0])
        vial_assignments = self.bin_vials(spots_true, vials = self.vials, bin_lines = self.bin_lines)[1]
        spots_true.loc[(spots_true.x >= self.bin_lines[0]) & (spots_true.x <= self.bin_lines[-1]),'vial'] = vial_assignments

        spots_true.loc[:,'color'] = spots_true.vial.map(dict(zip(range(1,self.vials+1), self.color_list)))
        bins=40

        ## Setting plots for scatterplot overlay on a selected frame
        if self.debug: print('detector.parameter_testing: Subplot 1: Test frame')
        axes[1].set_title('Frame: '+str(self.check_frame))
        axes[1].imshow(self.clean_stack[self.check_frame], cmap = cm.Greys_r)
        axes[1].scatter(spots_false[(spots_false.frame==self.check_frame)].x,
                        spots_false[(spots_false.frame==self.check_frame)].y, 
                        color = 'b',marker ='+',alpha = .5)
        a = axes[1].scatter(spots_true[spots_true.frame==self.check_frame].x,
                            spots_true[spots_true.frame==self.check_frame].y, 
                            c = spots_true[spots_true.frame==self.check_frame].vial,
                            cmap = self.vial_color_map,
                            marker ='o',alpha = .8)
        a.set_facecolor('none')
        axes[1].vlines(self.bin_lines,0,self.df_big.y.max(),color='w')
        axes[1].set_xlim(0,self.w)
        axes[1].set_ylim(self.h,0)
        
        ##########
        ## Fly counts
#         axes[5].plot(spots_true.groupby('frame').frame.mean(),
#                      spots_true.groupby('frame').frame.count(), 
#                      label = 'Fly count', color = 'g',alpha = .5)
#                      
#         axes[5].hlines(np.median(spots_true.groupby('frame').frame.count()),
#                        self.df_big.frame.min(),self.df_big.frame.max(),
#                        linestyle = '--',alpha = .5, 
#                        color = 'gray', label = 'Median no. flies')
#         axes[5].set(title = 'Flies per frame',
#                         ylabel='Flies detected',
#                         xlabel='Frame') 
#         axes[5].legend(frameon=False, fontsize = 'small')

        ##########
        df = self.df_filtered.sort_values(by='frame')
        for V in range(1,self.vials + 1):
            color = self.color_list[V-1]
            _df = df[df.vial == V]

            ## Local linear regression
            begin = self.local_linear_regression(_df).iloc[0].first_frame.astype(int)
            end = begin + self.window        
                        
            ## Plotting all points
            axes[5].plot(_df.groupby('frame').frame.unique(),
                _df.groupby('frame').y.count(),alpha = .3, color = color,label='') 
            
            ## Plotting most linear points
            _df = _df[(_df.frame >= begin) & (_df.frame <= end)]
            axes[5].plot(_df.groupby('frame').frame.unique(),
                _df.groupby('frame').frame.count() ,color = color, alpha = .5)
    
            axes[5].hlines(np.median(_df.groupby('frame').frame.count()),
                       df.frame.min(),df.frame.max(),
                       linestyle = '--',alpha = .7, 
                       color = color)


        # Deciding number of columns for legend
        if self.vials > 10: ncol = 3
        elif self.vials > 5: ncol = 2
        else: ncol=1
        
        ## Setting labels
        label_y,label_x = '(pixels)','Frames'
        if self.convert_to_cm_sec: 
            label_x,label_y = 'Seconds','(cm)'
        labels = ['Flies detected per frame','Flies detected','Frame']
        axes[5].set(title = labels[0], ylabel=labels[1],xlabel=labels[2]) 
        axes[5].set_ylim(ymin = 0,ymax = np.max(_df.groupby('frame').frame.count())*1.2)
#         axes[5].legend(frameon=False, fontsize='x-small', ncol=ncol)

        custom_lines = [Line2D([0], [0], color='k', linestyle = '--', alpha = .9),
                        Line2D([0], [0], color='k', linestyle = '-', alpha = .5)]
        custom_labels = ['Median', 'All frames']
        axes[5].legend(custom_lines, custom_labels,frameon=False, fontsize='x-small', ncol=ncol)

        #############



        ## Mass histogram
        axes[3].set_title('Mass Distribution')
        axes[3].hist(self.df_big.mass,bins = bins)
        y_max = np.histogram(self.df_big.mass,bins=bins)[0].max()
        axes[3].vlines(self.minmass,0,y_max)

        ## Signal histogram
        axes[4].set_title('Signal Distribution')
        axes[4].hist(self.df_big.signal,bins = bins)
        y_max = np.histogram(self.df_big.signal,bins=bins)[0].max()
        axes[4].vlines(self.threshold,0,y_max)

        ## Calculating local linear regression
#         self.step_5()
        
        ## Setting plots for local linear regression
        df = self.df_filtered.sort_values(by='frame')

        ## Converting to cm per sec if specified
        convert_x,convert_y = 1,1
        if self.convert_to_cm_sec:
            convert_x,convert_y = self.frame_rate,self.pixel_to_cm

        ## LocLin plot for each vial
        for V in range(1,self.vials + 1):
            label = 'Vial '+str(V)
            color = self.color_list[V-1]
            _df = df[df.vial == V]

            ## Local linear regression
            begin = self.local_linear_regression(_df).iloc[0].first_frame.astype(int)
            end = begin + self.window        
            
            ## Plotting all points
            axes[2].plot(_df.groupby('frame').frame.mean() / convert_x,
               _df.groupby('frame').y.mean() / convert_y,alpha = .35, color = color,label='') 
            
            ## Plotting most linear points
            _df = _df[(_df.frame >= begin) & (_df.frame <= end)]
            axes[2].plot(_df.groupby('frame').frame.mean() / convert_x,
           _df.groupby('frame').y.mean() / convert_y,color = color, label = label)

        # Deciding number of columns for legend
        if self.vials > 10: ncol = 3
        elif self.vials > 5: ncol = 2
        else: ncol=1
        
        ## Setting labels
        label_y,label_x = '(pixels)','Frames'
        if self.convert_to_cm_sec: 
            label_x,label_y = 'Seconds','(cm)'
        labels = ['Mean vertical position over time','Mean y-position %s' % label_y,label_x]
        axes[2].set(title = labels[0], ylabel=labels[1],xlabel=labels[2]) 
        axes[2].legend(frameon=False, fontsize='x-small', ncol=ncol)

        return