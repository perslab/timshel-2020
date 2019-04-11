############### SYNOPSIS ###################
# AIM: Make GIF/MP4 annitation of cell-type prioritization main figure


# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

library(gganimate)

setwd(here("src/publication"))


# ======================================================================= #
# ============================ TODO =============================== #
# ======================================================================= #

# mousebrain animation needs to be adjusted in p.main.margin. See Tabula Muris for how to solve these issues

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
dataset_prefix <- "tabula_muris"
# dataset_prefix <- "mousebrain_all"

# ======================================================================= #
# ============================ DEPENDENCIES =============================== #
# ======================================================================= #

# - p.main.margin: plot obtained from fig-celltype_priori_{mb/tm}.R . This is the main plot without any heatmaps etc.

### Get p.main.margin
if (dataset_prefix == "mousebrain_all") {
  source("fig-celltype_priori_mb.R")
  transition_state_variable <- sym("TaxonomyRank4_reduced1")
} else if (dataset_prefix == "tabula_muris") {
  source("fig-celltype_priori_tm.R")
  transition_state_variable <- sym("tissue")
}

# ======================================================================= #
# ============================ GGANIMATE =============================== #
# ======================================================================= #

### Create animation
p.ani <- p.main.margin + 
  transition_states(!!transition_state_variable,
                    transition_length=1, # The relative length of the transition
                    state_length=3, # The relative length of the pause at the states.
                    wrap = TRUE # You can set it to false if you don’t want the last state to transition back to the first state (default == TRUE).
                    # ^ becasue of a current BUG in gganimate wrap must be set to TRUE. Else you get: Error in data.frame(..., check.names = FALSE) :  arguments imply differing number of rows: 229, 100
                    # ^ BUG REF: https://github.com/thomasp85/gganimate/issues/301
  ) +
  shadow_mark() # + keep points visable
# enter_fade() # play with enter

### Add 'frame data' (mostly for debugging)
p.ani <- p.ani + labs(title = "Closest State: {closest_state}", 
                  subtitle = "Frame {frame} of {nframes}")

### Errros
# p.main.margin + transition_reveal(annotation) # Error: along data must either be integer, numeric, POSIXct, Date, difftime, orhms
# p.main.margin + transition_states(annotation, wrap = FALSE) + shadow_mark() # Error in data.frame(..., check.names = FALSE) :  arguments imply differing number of rows: 229, 100

flag_render_video <- TRUE
# flag_render_video <- FALSE

### RENDER video
# DURATION of this video: ~10 seconds works very well.
p.ani.rend <- animate(p.ani,
                      nframes=50, # number of frames MUST be more than number of states # default=100 | The length and framerate is decided on render time and can be any two combination of nframes, fps, and duration. 
                      # number of frames for transition_states()?: n_groups* (1 states + 1 transition) = 2n (I think?)
                      # fps=10, # default=10
                      duration=8, # seconds | default is unset
                      # ^ duration=nframes/fps # The length of the animation in seconds
                      detail=1, # default=1 | detail will get multiplied to nframes and the resulting number of frames will get calculated, but only nframes evenly spaced frames are rendered.
                      rewind=F,
                      # renderer=gifski_renderer(loop=F)
                      renderer=ifelse(flag_render_video, 
                                      av_renderer(), # do not use ffmpeg_renderer() because it ignores the your fps argument and only use fps=25, 
                                      gifski_renderer(loop=F) # no loop. PowerPoint2016 apparently ignores the looping metadata for any of my Gifs
                      ),
                      # end_pause=999, # Number of times to repeat the last frame in the animation | TRICKY TO GET TO WORK because it cannot exceed number of frames. It was an attempt to solving OSX PowerPoint looping issue with loop. 
                      # device options:
                      device="png", # default
                      width=9, 
                      height=8, 
                      units="in", 
                      res=150 # resolution in DPI
                      # width = 1280, height = 720, res = 104,
)


file.out.animation <- sprintf("figs/fig_animation.%s.%s", 
                              dataset_prefix,
                              ifelse(flag_render_video, "mp4", "gif"))
anim_save(filename=file.out.animation, animation=p.ani.rend)

# ======================= WIKI ======================= #

### pdf:
# width, height: the width and height of the graphics region in **inches**. The default values are 7.
### png: 
# units: the units in which height and width are given. Can be px (pixels, the *default*), in (inches), cm or mm.
# res: the nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer. Also used for units other than the default, and to convert points to pixels.
# animate(p.ani, detail=XXX)

### Gifsicle for editing gifs
# ^ See Gifsicle for command-line tool for GIF: https://www.lcdf.org/gifsicle/
# ^ gifsicle loop.gif --no-loopcount > no_loop.gif | REF: https://davidwalsh.name/prevent-gif-loop