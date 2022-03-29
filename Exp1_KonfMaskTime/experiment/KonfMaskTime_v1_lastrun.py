#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2020.2.10),
    on Juli 26, 2021, at 12:21
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2020.2.10'
expName = 'KonfMaskTime_v1'  # from the Builder filename that created this script
expInfo = {'gender': 'sag ich nicht', 'age': 'sag ich nicht', 'session': '', 'participant': ''}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s_%s' %(expInfo['participant'],expInfo['session'],  expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\PPA859\\Documents\\SeqSamplingConfidence\\KonfMaskTime\\KonfMaskTime_v1_lastrun.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1080], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True)
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "Instruktion"
InstruktionClock = core.Clock()
instruktion = visual.TextStim(win=win, name='instruktion',
    text='Herzlich Willkommen zu unserem Experiment!\nIn jedem Durchgang entscheiden Sie bitte, ob die Streifen des Gitters waagrecht oder senkrecht sind.\nZusätzlich geben Sie bitte an, wie sicher Sie sich dabei sind. \nWenn Sie keine weiteren Fragen mehr haben, drücken Sie bitte mit Ihrem Zeigefinger die Taste am Joystick, um zu beginnen.',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-1.0);

# Initialize components for Routine "trial"
trialClock = core.Clock()
from psychopy.hardware import joystick
joystick.backend='pyglet'
joy = joystick.Joystick(0)
newPosition = [0,0]
buttonNotYetPressed = True

ISI = clock.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI')
Fixationskreuz = visual.TextStim(win=win, name='Fixationskreuz',
    text='+',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-5.0);
TargetStimulus = visual.ImageStim(
    win=win,
    name='TargetStimulus', units='deg', 
    image='LowerContrastTarget.png', mask=None,
    ori=1.0, pos=[0, 0], size=[3, 3],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-6.0)
Maske = visual.ImageStim(
    win=win,
    name='Maske', units='deg', 
    image='sin', mask=None,
    ori=0, pos=[0, 0], size=[4, 4],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-7.0)
SkalaSenkrecht = visual.Rect(
    win=win, name='SkalaSenkrecht',units='deg', 
    width=[20, .1][0], height=[20, .1][1],
    ori=0, pos=[0, 10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-8.0, interpolate=True)
SkalaWaagrecht = visual.Rect(
    win=win, name='SkalaWaagrecht',units='deg', 
    width=[20, .1][0], height=[20, .1][1],
    ori=0, pos=[0, -10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-9.0, interpolate=True)
Marker = visual.Rect(
    win=win, name='Marker',units='deg', 
    width=[0.7, 0.7][0], height=[0.7, 0.7][1],
    ori=45, pos=[0,0],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-10.0, interpolate=True)
StrichObenLinks = visual.Rect(
    win=win, name='StrichObenLinks',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[-10, 10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-11.0, interpolate=True)
StrichObenRechts = visual.Rect(
    win=win, name='StrichObenRechts',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[10, 10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-12.0, interpolate=True)
StrichUntenLinks = visual.Rect(
    win=win, name='StrichUntenLinks',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[-10, -10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-13.0, interpolate=True)
StrichUntenRechts = visual.Rect(
    win=win, name='StrichUntenRechts',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[10, -10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-14.0, interpolate=True)
Senkrecht = visual.TextStim(win=win, name='Senkrecht',
    text='Senkrecht',
    font='Arial',
    units='deg', pos=[0, 9], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-15.0);
Waagrecht = visual.TextStim(win=win, name='Waagrecht',
    text='Waagrecht',
    font='Arial',
    units='deg', pos=[0, -9], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-16.0);
SicherOben = visual.TextStim(win=win, name='SicherOben',
    text='Sicher\nsenkrecht',
    font='Arial',
    units='deg', pos=[13, 10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-17.0);
SicherUnten = visual.TextStim(win=win, name='SicherUnten',
    text='Sicher\nwaagrecht',
    font='Arial',
    units='deg', pos=[13, -10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-18.0);
UnsicherOben = visual.TextStim(win=win, name='UnsicherOben',
    text='Unsicher\nsenkrecht',
    font='Arial',
    units='deg', pos=[-13, 10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-19.0);
UnsicherUnten = visual.TextStim(win=win, name='UnsicherUnten',
    text='Unsicher\nwaagrecht',
    font='Arial',
    units='deg', pos=[-13, -10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-20.0);

# Initialize components for Routine "Fehlerfeedback"
FehlerfeedbackClock = core.Clock()
Fehler = visual.TextStim(win=win, name='Fehler',
    text='Fehler!',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Initialize components for Routine "Blockende_Training"
Blockende_TrainingClock = core.Clock()
Durchgangsende = visual.TextStim(win=win, name='Durchgangsende',
    text='default text',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-2.0);

# Initialize components for Routine "trial"
trialClock = core.Clock()
from psychopy.hardware import joystick
joystick.backend='pyglet'
joy = joystick.Joystick(0)
newPosition = [0,0]
buttonNotYetPressed = True

ISI = clock.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI')
Fixationskreuz = visual.TextStim(win=win, name='Fixationskreuz',
    text='+',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-5.0);
TargetStimulus = visual.ImageStim(
    win=win,
    name='TargetStimulus', units='deg', 
    image='LowerContrastTarget.png', mask=None,
    ori=1.0, pos=[0, 0], size=[3, 3],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-6.0)
Maske = visual.ImageStim(
    win=win,
    name='Maske', units='deg', 
    image='sin', mask=None,
    ori=0, pos=[0, 0], size=[4, 4],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-7.0)
SkalaSenkrecht = visual.Rect(
    win=win, name='SkalaSenkrecht',units='deg', 
    width=[20, .1][0], height=[20, .1][1],
    ori=0, pos=[0, 10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-8.0, interpolate=True)
SkalaWaagrecht = visual.Rect(
    win=win, name='SkalaWaagrecht',units='deg', 
    width=[20, .1][0], height=[20, .1][1],
    ori=0, pos=[0, -10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-9.0, interpolate=True)
Marker = visual.Rect(
    win=win, name='Marker',units='deg', 
    width=[0.7, 0.7][0], height=[0.7, 0.7][1],
    ori=45, pos=[0,0],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-10.0, interpolate=True)
StrichObenLinks = visual.Rect(
    win=win, name='StrichObenLinks',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[-10, 10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-11.0, interpolate=True)
StrichObenRechts = visual.Rect(
    win=win, name='StrichObenRechts',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[10, 10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-12.0, interpolate=True)
StrichUntenLinks = visual.Rect(
    win=win, name='StrichUntenLinks',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[-10, -10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-13.0, interpolate=True)
StrichUntenRechts = visual.Rect(
    win=win, name='StrichUntenRechts',units='deg', 
    width=[0.1, 1][0], height=[0.1, 1][1],
    ori=0, pos=[10, -10],
    lineWidth=1, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-14.0, interpolate=True)
Senkrecht = visual.TextStim(win=win, name='Senkrecht',
    text='Senkrecht',
    font='Arial',
    units='deg', pos=[0, 9], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-15.0);
Waagrecht = visual.TextStim(win=win, name='Waagrecht',
    text='Waagrecht',
    font='Arial',
    units='deg', pos=[0, -9], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-16.0);
SicherOben = visual.TextStim(win=win, name='SicherOben',
    text='Sicher\nsenkrecht',
    font='Arial',
    units='deg', pos=[13, 10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-17.0);
SicherUnten = visual.TextStim(win=win, name='SicherUnten',
    text='Sicher\nwaagrecht',
    font='Arial',
    units='deg', pos=[13, -10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-18.0);
UnsicherOben = visual.TextStim(win=win, name='UnsicherOben',
    text='Unsicher\nsenkrecht',
    font='Arial',
    units='deg', pos=[-13, 10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-19.0);
UnsicherUnten = visual.TextStim(win=win, name='UnsicherUnten',
    text='Unsicher\nwaagrecht',
    font='Arial',
    units='deg', pos=[-13, -10], height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-20.0);

# Initialize components for Routine "Fehlerfeedback"
FehlerfeedbackClock = core.Clock()
Fehler = visual.TextStim(win=win, name='Fehler',
    text='Fehler!',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Initialize components for Routine "Blockende"
BlockendeClock = core.Clock()
Blockende_Text = visual.TextStim(win=win, name='Blockende_Text',
    text='Ende des Durchgangs\nBitte druecken Sie die Joysticktaste, um weiter zu machen',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=-2.0);

# Initialize components for Routine "Ende"
EndeClock = core.Clock()
Danke = visual.TextStim(win=win, name='Danke',
    text='Geschafft! :)\nVielen Dank für Ihre Teilnahme.\nMelden Sie sich bitte beim Versuchsleiter.',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "Instruktion"-------
continueRoutine = True
# update component parameters for each repeat
# keep track of which components have finished
InstruktionComponents = [instruktion]
for thisComponent in InstruktionComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
InstruktionClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "Instruktion"-------
while continueRoutine:
    # get current time
    t = InstruktionClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=InstruktionClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    buttonIsPressed = joy.getButton(0)
    
    if buttonIsPressed and buttonNotYetPressed:
     
      buttonNotYetPressed = False
    
    if ((not buttonIsPressed) and (not buttonNotYetPressed)):
        continueRoutine = False 
    
    # *instruktion* updates
    if instruktion.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        instruktion.frameNStart = frameN  # exact frame index
        instruktion.tStart = t  # local t and not account for scr refresh
        instruktion.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(instruktion, 'tStartRefresh')  # time at next scr refresh
        instruktion.setAutoDraw(True)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in InstruktionComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "Instruktion"-------
for thisComponent in InstruktionComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('instruktion.started', instruktion.tStartRefresh)
thisExp.addData('instruktion.stopped', instruktion.tStopRefresh)
# the Routine "Instruktion" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
trials_Training = data.TrialHandler(nReps=12, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('BedingungenMasking.xlsx'),
    seed=None, name='trials_Training')
thisExp.addLoop(trials_Training)  # add the loop to the experiment
thisTrials_Training = trials_Training.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrials_Training.rgb)
if thisTrials_Training != None:
    for paramName in thisTrials_Training:
        exec('{} = thisTrials_Training[paramName]'.format(paramName))

for thisTrials_Training in trials_Training:
    currentLoop = trials_Training
    # abbreviate parameter names if possible (e.g. rgb = thisTrials_Training.rgb)
    if thisTrials_Training != None:
        for paramName in thisTrials_Training:
            exec('{} = thisTrials_Training[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "trial"-------
    continueRoutine = True
    # update component parameters for each repeat
    ori = np.random.choice([0,90])
    thisExp.addData('Orientrierung', ori)
    
    if ori == 90:
        correctAns = 'waagrecht'
    else:
        correctAns = 'senkrecht'
    
    thisExp.addData('expectedAnswer', correctAns)
    buttonNotYetPressed = True
    
    TargetStimulus.setOri(ori)
    Maske.setImage('checkerBoardSmall.png')
    # keep track of which components have finished
    trialComponents = [ISI, Fixationskreuz, TargetStimulus, Maske, SkalaSenkrecht, SkalaWaagrecht, Marker, StrichObenLinks, StrichObenRechts, StrichUntenLinks, StrichUntenRechts, Senkrecht, Waagrecht, SicherOben, SicherUnten, UnsicherOben, UnsicherUnten]
    for thisComponent in trialComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "trial"-------
    while continueRoutine:
        # get current time
        t = trialClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=trialClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        xPos = joy.getX()
        joyyPos = joy.getY()
        
        if (joyyPos <= -0.70):
          #showMarker = True
          yPos = 10
          ans = 'senkrecht'
        
        elif (joyyPos >= 0.70):
          #showMarker = True
          yPos = - 10
          ans = 'waagrecht'
        
        else: 
         #showMarker = False
         yPos = 20000
        
        newPosition = [xPos * 10, yPos]
        buttonIsPressed = joy.getButton(0)
        
        if buttonIsPressed and buttonNotYetPressed and (trialClock.getTime() >= (1.5+SOA/1000)) and (yPos != 20000):
          currentTime = trialClock.getTime()
          thisExp.addData('Ans', ans)
          thisExp.addData('Rating', xPos)
          thisExp.addData('RatingRT', currentTime)
          if (correctAns == ans):
            correct = 1
          else: 
            correct = 0
          thisExp.addData('ResponseCorr', correct)
          buttonNotYetPressed = False
        
        if ((not buttonIsPressed) and (not buttonNotYetPressed)):
            continueRoutine = False 
        
        # *Fixationskreuz* updates
        if Fixationskreuz.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            Fixationskreuz.frameNStart = frameN  # exact frame index
            Fixationskreuz.tStart = t  # local t and not account for scr refresh
            Fixationskreuz.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Fixationskreuz, 'tStartRefresh')  # time at next scr refresh
            Fixationskreuz.setAutoDraw(True)
        if Fixationskreuz.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > Fixationskreuz.tStartRefresh + 1.0-frameTolerance:
                # keep track of stop time/frame for later
                Fixationskreuz.tStop = t  # not accounting for scr refresh
                Fixationskreuz.frameNStop = frameN  # exact frame index
                win.timeOnFlip(Fixationskreuz, 'tStopRefresh')  # time at next scr refresh
                Fixationskreuz.setAutoDraw(False)
        
        # *TargetStimulus* updates
        if TargetStimulus.status == NOT_STARTED and tThisFlip >= 1-frameTolerance:
            # keep track of start time/frame for later
            TargetStimulus.frameNStart = frameN  # exact frame index
            TargetStimulus.tStart = t  # local t and not account for scr refresh
            TargetStimulus.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(TargetStimulus, 'tStartRefresh')  # time at next scr refresh
            TargetStimulus.setAutoDraw(True)
        if TargetStimulus.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > TargetStimulus.tStartRefresh + SOA/1000-frameTolerance:
                # keep track of stop time/frame for later
                TargetStimulus.tStop = t  # not accounting for scr refresh
                TargetStimulus.frameNStop = frameN  # exact frame index
                win.timeOnFlip(TargetStimulus, 'tStopRefresh')  # time at next scr refresh
                TargetStimulus.setAutoDraw(False)
        
        # *Maske* updates
        if Maske.status == NOT_STARTED and tThisFlip >= 1+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            Maske.frameNStart = frameN  # exact frame index
            Maske.tStart = t  # local t and not account for scr refresh
            Maske.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Maske, 'tStartRefresh')  # time at next scr refresh
            Maske.setAutoDraw(True)
        if Maske.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > Maske.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                Maske.tStop = t  # not accounting for scr refresh
                Maske.frameNStop = frameN  # exact frame index
                win.timeOnFlip(Maske, 'tStopRefresh')  # time at next scr refresh
                Maske.setAutoDraw(False)
        
        # *SkalaSenkrecht* updates
        if SkalaSenkrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            SkalaSenkrecht.frameNStart = frameN  # exact frame index
            SkalaSenkrecht.tStart = t  # local t and not account for scr refresh
            SkalaSenkrecht.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(SkalaSenkrecht, 'tStartRefresh')  # time at next scr refresh
            SkalaSenkrecht.setAutoDraw(True)
        
        # *SkalaWaagrecht* updates
        if SkalaWaagrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            SkalaWaagrecht.frameNStart = frameN  # exact frame index
            SkalaWaagrecht.tStart = t  # local t and not account for scr refresh
            SkalaWaagrecht.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(SkalaWaagrecht, 'tStartRefresh')  # time at next scr refresh
            SkalaWaagrecht.setAutoDraw(True)
        
        # *Marker* updates
        if Marker.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            Marker.frameNStart = frameN  # exact frame index
            Marker.tStart = t  # local t and not account for scr refresh
            Marker.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Marker, 'tStartRefresh')  # time at next scr refresh
            Marker.setAutoDraw(True)
        if Marker.status == STARTED:  # only update if drawing
            Marker.setPos(newPosition)
        
        # *StrichObenLinks* updates
        if StrichObenLinks.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            StrichObenLinks.frameNStart = frameN  # exact frame index
            StrichObenLinks.tStart = t  # local t and not account for scr refresh
            StrichObenLinks.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(StrichObenLinks, 'tStartRefresh')  # time at next scr refresh
            StrichObenLinks.setAutoDraw(True)
        
        # *StrichObenRechts* updates
        if StrichObenRechts.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            StrichObenRechts.frameNStart = frameN  # exact frame index
            StrichObenRechts.tStart = t  # local t and not account for scr refresh
            StrichObenRechts.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(StrichObenRechts, 'tStartRefresh')  # time at next scr refresh
            StrichObenRechts.setAutoDraw(True)
        
        # *StrichUntenLinks* updates
        if StrichUntenLinks.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            StrichUntenLinks.frameNStart = frameN  # exact frame index
            StrichUntenLinks.tStart = t  # local t and not account for scr refresh
            StrichUntenLinks.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(StrichUntenLinks, 'tStartRefresh')  # time at next scr refresh
            StrichUntenLinks.setAutoDraw(True)
        
        # *StrichUntenRechts* updates
        if StrichUntenRechts.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            StrichUntenRechts.frameNStart = frameN  # exact frame index
            StrichUntenRechts.tStart = t  # local t and not account for scr refresh
            StrichUntenRechts.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(StrichUntenRechts, 'tStartRefresh')  # time at next scr refresh
            StrichUntenRechts.setAutoDraw(True)
        
        # *Senkrecht* updates
        if Senkrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            Senkrecht.frameNStart = frameN  # exact frame index
            Senkrecht.tStart = t  # local t and not account for scr refresh
            Senkrecht.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Senkrecht, 'tStartRefresh')  # time at next scr refresh
            Senkrecht.setAutoDraw(True)
        
        # *Waagrecht* updates
        if Waagrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            Waagrecht.frameNStart = frameN  # exact frame index
            Waagrecht.tStart = t  # local t and not account for scr refresh
            Waagrecht.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Waagrecht, 'tStartRefresh')  # time at next scr refresh
            Waagrecht.setAutoDraw(True)
        
        # *SicherOben* updates
        if SicherOben.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            SicherOben.frameNStart = frameN  # exact frame index
            SicherOben.tStart = t  # local t and not account for scr refresh
            SicherOben.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(SicherOben, 'tStartRefresh')  # time at next scr refresh
            SicherOben.setAutoDraw(True)
        
        # *SicherUnten* updates
        if SicherUnten.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            SicherUnten.frameNStart = frameN  # exact frame index
            SicherUnten.tStart = t  # local t and not account for scr refresh
            SicherUnten.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(SicherUnten, 'tStartRefresh')  # time at next scr refresh
            SicherUnten.setAutoDraw(True)
        
        # *UnsicherOben* updates
        if UnsicherOben.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            UnsicherOben.frameNStart = frameN  # exact frame index
            UnsicherOben.tStart = t  # local t and not account for scr refresh
            UnsicherOben.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(UnsicherOben, 'tStartRefresh')  # time at next scr refresh
            UnsicherOben.setAutoDraw(True)
        
        # *UnsicherUnten* updates
        if UnsicherUnten.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
            # keep track of start time/frame for later
            UnsicherUnten.frameNStart = frameN  # exact frame index
            UnsicherUnten.tStart = t  # local t and not account for scr refresh
            UnsicherUnten.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(UnsicherUnten, 'tStartRefresh')  # time at next scr refresh
            UnsicherUnten.setAutoDraw(True)
        # *ISI* period
        if ISI.status == NOT_STARTED and t >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            ISI.frameNStart = frameN  # exact frame index
            ISI.tStart = t  # local t and not account for scr refresh
            ISI.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(ISI, 'tStartRefresh')  # time at next scr refresh
            ISI.start(0.5)
        elif ISI.status == STARTED:  # one frame should pass before updating params and completing
            ISI.complete()  # finish the static period
            ISI.tStop = ISI.tStart + 0.5  # record stop time
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "trial"-------
    for thisComponent in trialComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    trials_Training.addData('ISI.started', ISI.tStart)
    trials_Training.addData('ISI.stopped', ISI.tStop)
    trials_Training.addData('Fixationskreuz.started', Fixationskreuz.tStartRefresh)
    trials_Training.addData('Fixationskreuz.stopped', Fixationskreuz.tStopRefresh)
    trials_Training.addData('TargetStimulus.started', TargetStimulus.tStartRefresh)
    trials_Training.addData('TargetStimulus.stopped', TargetStimulus.tStopRefresh)
    trials_Training.addData('Maske.started', Maske.tStartRefresh)
    trials_Training.addData('Maske.stopped', Maske.tStopRefresh)
    trials_Training.addData('SkalaSenkrecht.started', SkalaSenkrecht.tStartRefresh)
    trials_Training.addData('SkalaSenkrecht.stopped', SkalaSenkrecht.tStopRefresh)
    trials_Training.addData('SkalaWaagrecht.started', SkalaWaagrecht.tStartRefresh)
    trials_Training.addData('SkalaWaagrecht.stopped', SkalaWaagrecht.tStopRefresh)
    trials_Training.addData('Marker.started', Marker.tStartRefresh)
    trials_Training.addData('Marker.stopped', Marker.tStopRefresh)
    trials_Training.addData('StrichObenLinks.started', StrichObenLinks.tStartRefresh)
    trials_Training.addData('StrichObenLinks.stopped', StrichObenLinks.tStopRefresh)
    trials_Training.addData('StrichObenRechts.started', StrichObenRechts.tStartRefresh)
    trials_Training.addData('StrichObenRechts.stopped', StrichObenRechts.tStopRefresh)
    trials_Training.addData('StrichUntenLinks.started', StrichUntenLinks.tStartRefresh)
    trials_Training.addData('StrichUntenLinks.stopped', StrichUntenLinks.tStopRefresh)
    trials_Training.addData('StrichUntenRechts.started', StrichUntenRechts.tStartRefresh)
    trials_Training.addData('StrichUntenRechts.stopped', StrichUntenRechts.tStopRefresh)
    trials_Training.addData('Senkrecht.started', Senkrecht.tStartRefresh)
    trials_Training.addData('Senkrecht.stopped', Senkrecht.tStopRefresh)
    trials_Training.addData('Waagrecht.started', Waagrecht.tStartRefresh)
    trials_Training.addData('Waagrecht.stopped', Waagrecht.tStopRefresh)
    trials_Training.addData('SicherOben.started', SicherOben.tStartRefresh)
    trials_Training.addData('SicherOben.stopped', SicherOben.tStopRefresh)
    trials_Training.addData('SicherUnten.started', SicherUnten.tStartRefresh)
    trials_Training.addData('SicherUnten.stopped', SicherUnten.tStopRefresh)
    trials_Training.addData('UnsicherOben.started', UnsicherOben.tStartRefresh)
    trials_Training.addData('UnsicherOben.stopped', UnsicherOben.tStopRefresh)
    trials_Training.addData('UnsicherUnten.started', UnsicherUnten.tStartRefresh)
    trials_Training.addData('UnsicherUnten.stopped', UnsicherUnten.tStopRefresh)
    # the Routine "trial" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    Feedback_Training = data.TrialHandler(nReps=1 - correct, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='Feedback_Training')
    thisExp.addLoop(Feedback_Training)  # add the loop to the experiment
    thisFeedback_Training = Feedback_Training.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisFeedback_Training.rgb)
    if thisFeedback_Training != None:
        for paramName in thisFeedback_Training:
            exec('{} = thisFeedback_Training[paramName]'.format(paramName))
    
    for thisFeedback_Training in Feedback_Training:
        currentLoop = Feedback_Training
        # abbreviate parameter names if possible (e.g. rgb = thisFeedback_Training.rgb)
        if thisFeedback_Training != None:
            for paramName in thisFeedback_Training:
                exec('{} = thisFeedback_Training[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "Fehlerfeedback"-------
        continueRoutine = True
        routineTimer.add(1.000000)
        # update component parameters for each repeat
        # keep track of which components have finished
        FehlerfeedbackComponents = [Fehler]
        for thisComponent in FehlerfeedbackComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        FehlerfeedbackClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "Fehlerfeedback"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = FehlerfeedbackClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=FehlerfeedbackClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *Fehler* updates
            if Fehler.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                Fehler.frameNStart = frameN  # exact frame index
                Fehler.tStart = t  # local t and not account for scr refresh
                Fehler.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Fehler, 'tStartRefresh')  # time at next scr refresh
                Fehler.setAutoDraw(True)
            if Fehler.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Fehler.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    Fehler.tStop = t  # not accounting for scr refresh
                    Fehler.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(Fehler, 'tStopRefresh')  # time at next scr refresh
                    Fehler.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in FehlerfeedbackComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "Fehlerfeedback"-------
        for thisComponent in FehlerfeedbackComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        Feedback_Training.addData('Fehler.started', Fehler.tStartRefresh)
        Feedback_Training.addData('Fehler.stopped', Fehler.tStopRefresh)
    # completed 1 - correct repeats of 'Feedback_Training'
    
    thisExp.nextEntry()
    
# completed 12 repeats of 'trials_Training'


# ------Prepare to start Routine "Blockende_Training"-------
continueRoutine = True
# update component parameters for each repeat
buttonNotYetPressed = True

Durchgangsende.setText('Ende des Durchgangs\nBitte druecken Sie die Joysticktaste, um weiter zu machen')
# keep track of which components have finished
Blockende_TrainingComponents = [Durchgangsende]
for thisComponent in Blockende_TrainingComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
Blockende_TrainingClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "Blockende_Training"-------
while continueRoutine:
    # get current time
    t = Blockende_TrainingClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=Blockende_TrainingClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    buttonIsPressed = joy.getButton(0)
    
    if buttonIsPressed and buttonNotYetPressed:
     
      buttonNotYetPressed = False
    
    if ((not buttonIsPressed) and (not buttonNotYetPressed)):
        continueRoutine = False 
    
    # *Durchgangsende* updates
    if Durchgangsende.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        Durchgangsende.frameNStart = frameN  # exact frame index
        Durchgangsende.tStart = t  # local t and not account for scr refresh
        Durchgangsende.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(Durchgangsende, 'tStartRefresh')  # time at next scr refresh
        Durchgangsende.setAutoDraw(True)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in Blockende_TrainingComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "Blockende_Training"-------
for thisComponent in Blockende_TrainingComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('Durchgangsende.started', Durchgangsende.tStartRefresh)
thisExp.addData('Durchgangsende.stopped', Durchgangsende.tStopRefresh)
# the Routine "Blockende_Training" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
blocks = data.TrialHandler(nReps=9, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='blocks')
thisExp.addLoop(blocks)  # add the loop to the experiment
thisBlock = blocks.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
if thisBlock != None:
    for paramName in thisBlock:
        exec('{} = thisBlock[paramName]'.format(paramName))

for thisBlock in blocks:
    currentLoop = blocks
    # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
    if thisBlock != None:
        for paramName in thisBlock:
            exec('{} = thisBlock[paramName]'.format(paramName))
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=12, method='fullRandom', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions('BedingungenMasking.xlsx'),
        seed=None, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            exec('{} = thisTrial[paramName]'.format(paramName))
    
    for thisTrial in trials:
        currentLoop = trials
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                exec('{} = thisTrial[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        # update component parameters for each repeat
        ori = np.random.choice([0,90])
        thisExp.addData('Orientrierung', ori)
        
        if ori == 90:
            correctAns = 'waagrecht'
        else:
            correctAns = 'senkrecht'
        
        thisExp.addData('expectedAnswer', correctAns)
        buttonNotYetPressed = True
        
        TargetStimulus.setOri(ori)
        Maske.setImage('checkerBoardSmall.png')
        # keep track of which components have finished
        trialComponents = [ISI, Fixationskreuz, TargetStimulus, Maske, SkalaSenkrecht, SkalaWaagrecht, Marker, StrichObenLinks, StrichObenRechts, StrichUntenLinks, StrichUntenRechts, Senkrecht, Waagrecht, SicherOben, SicherUnten, UnsicherOben, UnsicherUnten]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            xPos = joy.getX()
            joyyPos = joy.getY()
            
            if (joyyPos <= -0.70):
              #showMarker = True
              yPos = 10
              ans = 'senkrecht'
            
            elif (joyyPos >= 0.70):
              #showMarker = True
              yPos = - 10
              ans = 'waagrecht'
            
            else: 
             #showMarker = False
             yPos = 20000
            
            newPosition = [xPos * 10, yPos]
            buttonIsPressed = joy.getButton(0)
            
            if buttonIsPressed and buttonNotYetPressed and (trialClock.getTime() >= (1.5+SOA/1000)) and (yPos != 20000):
              currentTime = trialClock.getTime()
              thisExp.addData('Ans', ans)
              thisExp.addData('Rating', xPos)
              thisExp.addData('RatingRT', currentTime)
              if (correctAns == ans):
                correct = 1
              else: 
                correct = 0
              thisExp.addData('ResponseCorr', correct)
              buttonNotYetPressed = False
            
            if ((not buttonIsPressed) and (not buttonNotYetPressed)):
                continueRoutine = False 
            
            # *Fixationskreuz* updates
            if Fixationskreuz.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                Fixationskreuz.frameNStart = frameN  # exact frame index
                Fixationskreuz.tStart = t  # local t and not account for scr refresh
                Fixationskreuz.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Fixationskreuz, 'tStartRefresh')  # time at next scr refresh
                Fixationskreuz.setAutoDraw(True)
            if Fixationskreuz.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Fixationskreuz.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    Fixationskreuz.tStop = t  # not accounting for scr refresh
                    Fixationskreuz.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(Fixationskreuz, 'tStopRefresh')  # time at next scr refresh
                    Fixationskreuz.setAutoDraw(False)
            
            # *TargetStimulus* updates
            if TargetStimulus.status == NOT_STARTED and tThisFlip >= 1-frameTolerance:
                # keep track of start time/frame for later
                TargetStimulus.frameNStart = frameN  # exact frame index
                TargetStimulus.tStart = t  # local t and not account for scr refresh
                TargetStimulus.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(TargetStimulus, 'tStartRefresh')  # time at next scr refresh
                TargetStimulus.setAutoDraw(True)
            if TargetStimulus.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > TargetStimulus.tStartRefresh + SOA/1000-frameTolerance:
                    # keep track of stop time/frame for later
                    TargetStimulus.tStop = t  # not accounting for scr refresh
                    TargetStimulus.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(TargetStimulus, 'tStopRefresh')  # time at next scr refresh
                    TargetStimulus.setAutoDraw(False)
            
            # *Maske* updates
            if Maske.status == NOT_STARTED and tThisFlip >= 1+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                Maske.frameNStart = frameN  # exact frame index
                Maske.tStart = t  # local t and not account for scr refresh
                Maske.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Maske, 'tStartRefresh')  # time at next scr refresh
                Maske.setAutoDraw(True)
            if Maske.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Maske.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    Maske.tStop = t  # not accounting for scr refresh
                    Maske.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(Maske, 'tStopRefresh')  # time at next scr refresh
                    Maske.setAutoDraw(False)
            
            # *SkalaSenkrecht* updates
            if SkalaSenkrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                SkalaSenkrecht.frameNStart = frameN  # exact frame index
                SkalaSenkrecht.tStart = t  # local t and not account for scr refresh
                SkalaSenkrecht.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(SkalaSenkrecht, 'tStartRefresh')  # time at next scr refresh
                SkalaSenkrecht.setAutoDraw(True)
            
            # *SkalaWaagrecht* updates
            if SkalaWaagrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                SkalaWaagrecht.frameNStart = frameN  # exact frame index
                SkalaWaagrecht.tStart = t  # local t and not account for scr refresh
                SkalaWaagrecht.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(SkalaWaagrecht, 'tStartRefresh')  # time at next scr refresh
                SkalaWaagrecht.setAutoDraw(True)
            
            # *Marker* updates
            if Marker.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                Marker.frameNStart = frameN  # exact frame index
                Marker.tStart = t  # local t and not account for scr refresh
                Marker.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Marker, 'tStartRefresh')  # time at next scr refresh
                Marker.setAutoDraw(True)
            if Marker.status == STARTED:  # only update if drawing
                Marker.setPos(newPosition)
            
            # *StrichObenLinks* updates
            if StrichObenLinks.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                StrichObenLinks.frameNStart = frameN  # exact frame index
                StrichObenLinks.tStart = t  # local t and not account for scr refresh
                StrichObenLinks.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(StrichObenLinks, 'tStartRefresh')  # time at next scr refresh
                StrichObenLinks.setAutoDraw(True)
            
            # *StrichObenRechts* updates
            if StrichObenRechts.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                StrichObenRechts.frameNStart = frameN  # exact frame index
                StrichObenRechts.tStart = t  # local t and not account for scr refresh
                StrichObenRechts.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(StrichObenRechts, 'tStartRefresh')  # time at next scr refresh
                StrichObenRechts.setAutoDraw(True)
            
            # *StrichUntenLinks* updates
            if StrichUntenLinks.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                StrichUntenLinks.frameNStart = frameN  # exact frame index
                StrichUntenLinks.tStart = t  # local t and not account for scr refresh
                StrichUntenLinks.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(StrichUntenLinks, 'tStartRefresh')  # time at next scr refresh
                StrichUntenLinks.setAutoDraw(True)
            
            # *StrichUntenRechts* updates
            if StrichUntenRechts.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                StrichUntenRechts.frameNStart = frameN  # exact frame index
                StrichUntenRechts.tStart = t  # local t and not account for scr refresh
                StrichUntenRechts.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(StrichUntenRechts, 'tStartRefresh')  # time at next scr refresh
                StrichUntenRechts.setAutoDraw(True)
            
            # *Senkrecht* updates
            if Senkrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                Senkrecht.frameNStart = frameN  # exact frame index
                Senkrecht.tStart = t  # local t and not account for scr refresh
                Senkrecht.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Senkrecht, 'tStartRefresh')  # time at next scr refresh
                Senkrecht.setAutoDraw(True)
            
            # *Waagrecht* updates
            if Waagrecht.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                Waagrecht.frameNStart = frameN  # exact frame index
                Waagrecht.tStart = t  # local t and not account for scr refresh
                Waagrecht.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Waagrecht, 'tStartRefresh')  # time at next scr refresh
                Waagrecht.setAutoDraw(True)
            
            # *SicherOben* updates
            if SicherOben.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                SicherOben.frameNStart = frameN  # exact frame index
                SicherOben.tStart = t  # local t and not account for scr refresh
                SicherOben.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(SicherOben, 'tStartRefresh')  # time at next scr refresh
                SicherOben.setAutoDraw(True)
            
            # *SicherUnten* updates
            if SicherUnten.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                SicherUnten.frameNStart = frameN  # exact frame index
                SicherUnten.tStart = t  # local t and not account for scr refresh
                SicherUnten.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(SicherUnten, 'tStartRefresh')  # time at next scr refresh
                SicherUnten.setAutoDraw(True)
            
            # *UnsicherOben* updates
            if UnsicherOben.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                UnsicherOben.frameNStart = frameN  # exact frame index
                UnsicherOben.tStart = t  # local t and not account for scr refresh
                UnsicherOben.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(UnsicherOben, 'tStartRefresh')  # time at next scr refresh
                UnsicherOben.setAutoDraw(True)
            
            # *UnsicherUnten* updates
            if UnsicherUnten.status == NOT_STARTED and tThisFlip >= 1.5+SOA/1000-frameTolerance:
                # keep track of start time/frame for later
                UnsicherUnten.frameNStart = frameN  # exact frame index
                UnsicherUnten.tStart = t  # local t and not account for scr refresh
                UnsicherUnten.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(UnsicherUnten, 'tStartRefresh')  # time at next scr refresh
                UnsicherUnten.setAutoDraw(True)
            # *ISI* period
            if ISI.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                ISI.frameNStart = frameN  # exact frame index
                ISI.tStart = t  # local t and not account for scr refresh
                ISI.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(ISI, 'tStartRefresh')  # time at next scr refresh
                ISI.start(0.5)
            elif ISI.status == STARTED:  # one frame should pass before updating params and completing
                ISI.complete()  # finish the static period
                ISI.tStop = ISI.tStart + 0.5  # record stop time
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        trials.addData('ISI.started', ISI.tStart)
        trials.addData('ISI.stopped', ISI.tStop)
        trials.addData('Fixationskreuz.started', Fixationskreuz.tStartRefresh)
        trials.addData('Fixationskreuz.stopped', Fixationskreuz.tStopRefresh)
        trials.addData('TargetStimulus.started', TargetStimulus.tStartRefresh)
        trials.addData('TargetStimulus.stopped', TargetStimulus.tStopRefresh)
        trials.addData('Maske.started', Maske.tStartRefresh)
        trials.addData('Maske.stopped', Maske.tStopRefresh)
        trials.addData('SkalaSenkrecht.started', SkalaSenkrecht.tStartRefresh)
        trials.addData('SkalaSenkrecht.stopped', SkalaSenkrecht.tStopRefresh)
        trials.addData('SkalaWaagrecht.started', SkalaWaagrecht.tStartRefresh)
        trials.addData('SkalaWaagrecht.stopped', SkalaWaagrecht.tStopRefresh)
        trials.addData('Marker.started', Marker.tStartRefresh)
        trials.addData('Marker.stopped', Marker.tStopRefresh)
        trials.addData('StrichObenLinks.started', StrichObenLinks.tStartRefresh)
        trials.addData('StrichObenLinks.stopped', StrichObenLinks.tStopRefresh)
        trials.addData('StrichObenRechts.started', StrichObenRechts.tStartRefresh)
        trials.addData('StrichObenRechts.stopped', StrichObenRechts.tStopRefresh)
        trials.addData('StrichUntenLinks.started', StrichUntenLinks.tStartRefresh)
        trials.addData('StrichUntenLinks.stopped', StrichUntenLinks.tStopRefresh)
        trials.addData('StrichUntenRechts.started', StrichUntenRechts.tStartRefresh)
        trials.addData('StrichUntenRechts.stopped', StrichUntenRechts.tStopRefresh)
        trials.addData('Senkrecht.started', Senkrecht.tStartRefresh)
        trials.addData('Senkrecht.stopped', Senkrecht.tStopRefresh)
        trials.addData('Waagrecht.started', Waagrecht.tStartRefresh)
        trials.addData('Waagrecht.stopped', Waagrecht.tStopRefresh)
        trials.addData('SicherOben.started', SicherOben.tStartRefresh)
        trials.addData('SicherOben.stopped', SicherOben.tStopRefresh)
        trials.addData('SicherUnten.started', SicherUnten.tStartRefresh)
        trials.addData('SicherUnten.stopped', SicherUnten.tStopRefresh)
        trials.addData('UnsicherOben.started', UnsicherOben.tStartRefresh)
        trials.addData('UnsicherOben.stopped', UnsicherOben.tStopRefresh)
        trials.addData('UnsicherUnten.started', UnsicherUnten.tStartRefresh)
        trials.addData('UnsicherUnten.stopped', UnsicherUnten.tStopRefresh)
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        WarEsEinFehler = data.TrialHandler(nReps=1 - correct, method='random', 
            extraInfo=expInfo, originPath=-1,
            trialList=[None],
            seed=None, name='WarEsEinFehler')
        thisExp.addLoop(WarEsEinFehler)  # add the loop to the experiment
        thisWarEsEinFehler = WarEsEinFehler.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisWarEsEinFehler.rgb)
        if thisWarEsEinFehler != None:
            for paramName in thisWarEsEinFehler:
                exec('{} = thisWarEsEinFehler[paramName]'.format(paramName))
        
        for thisWarEsEinFehler in WarEsEinFehler:
            currentLoop = WarEsEinFehler
            # abbreviate parameter names if possible (e.g. rgb = thisWarEsEinFehler.rgb)
            if thisWarEsEinFehler != None:
                for paramName in thisWarEsEinFehler:
                    exec('{} = thisWarEsEinFehler[paramName]'.format(paramName))
            
            # ------Prepare to start Routine "Fehlerfeedback"-------
            continueRoutine = True
            routineTimer.add(1.000000)
            # update component parameters for each repeat
            # keep track of which components have finished
            FehlerfeedbackComponents = [Fehler]
            for thisComponent in FehlerfeedbackComponents:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            FehlerfeedbackClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
            frameN = -1
            
            # -------Run Routine "Fehlerfeedback"-------
            while continueRoutine and routineTimer.getTime() > 0:
                # get current time
                t = FehlerfeedbackClock.getTime()
                tThisFlip = win.getFutureFlipTime(clock=FehlerfeedbackClock)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                
                # *Fehler* updates
                if Fehler.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    Fehler.frameNStart = frameN  # exact frame index
                    Fehler.tStart = t  # local t and not account for scr refresh
                    Fehler.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(Fehler, 'tStartRefresh')  # time at next scr refresh
                    Fehler.setAutoDraw(True)
                if Fehler.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > Fehler.tStartRefresh + 1.0-frameTolerance:
                        # keep track of stop time/frame for later
                        Fehler.tStop = t  # not accounting for scr refresh
                        Fehler.frameNStop = frameN  # exact frame index
                        win.timeOnFlip(Fehler, 'tStopRefresh')  # time at next scr refresh
                        Fehler.setAutoDraw(False)
                
                # check for quit (typically the Esc key)
                if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                    core.quit()
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in FehlerfeedbackComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # -------Ending Routine "Fehlerfeedback"-------
            for thisComponent in FehlerfeedbackComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            WarEsEinFehler.addData('Fehler.started', Fehler.tStartRefresh)
            WarEsEinFehler.addData('Fehler.stopped', Fehler.tStopRefresh)
        # completed 1 - correct repeats of 'WarEsEinFehler'
        
        thisExp.nextEntry()
        
    # completed 12 repeats of 'trials'
    
    
    # ------Prepare to start Routine "Blockende"-------
    continueRoutine = True
    # update component parameters for each repeat
    buttonNotYetPressed = True
    
    # keep track of which components have finished
    BlockendeComponents = [Blockende_Text]
    for thisComponent in BlockendeComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    BlockendeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "Blockende"-------
    while continueRoutine:
        # get current time
        t = BlockendeClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=BlockendeClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        buttonIsPressed = joy.getButton(0)
        
        if buttonIsPressed and buttonNotYetPressed:
         
          buttonNotYetPressed = False
        
        if ((not buttonIsPressed) and (not buttonNotYetPressed)):
            continueRoutine = False 
        
        # *Blockende_Text* updates
        if Blockende_Text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            Blockende_Text.frameNStart = frameN  # exact frame index
            Blockende_Text.tStart = t  # local t and not account for scr refresh
            Blockende_Text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Blockende_Text, 'tStartRefresh')  # time at next scr refresh
            Blockende_Text.setAutoDraw(True)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in BlockendeComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "Blockende"-------
    for thisComponent in BlockendeComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    blocks.addData('Blockende_Text.started', Blockende_Text.tStartRefresh)
    blocks.addData('Blockende_Text.stopped', Blockende_Text.tStopRefresh)
    # the Routine "Blockende" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
# completed 9 repeats of 'blocks'


# ------Prepare to start Routine "Ende"-------
continueRoutine = True
routineTimer.add(5.000000)
# update component parameters for each repeat
# keep track of which components have finished
EndeComponents = [Danke]
for thisComponent in EndeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
EndeClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "Ende"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = EndeClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=EndeClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *Danke* updates
    if Danke.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        Danke.frameNStart = frameN  # exact frame index
        Danke.tStart = t  # local t and not account for scr refresh
        Danke.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(Danke, 'tStartRefresh')  # time at next scr refresh
        Danke.setAutoDraw(True)
    if Danke.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > Danke.tStartRefresh + 5-frameTolerance:
            # keep track of stop time/frame for later
            Danke.tStop = t  # not accounting for scr refresh
            Danke.frameNStop = frameN  # exact frame index
            win.timeOnFlip(Danke, 'tStopRefresh')  # time at next scr refresh
            Danke.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in EndeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "Ende"-------
for thisComponent in EndeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('Danke.started', Danke.tStartRefresh)
thisExp.addData('Danke.stopped', Danke.tStopRefresh)

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
