---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# Metronome
**Maintainer**: Simon Spannagel (<simon.spannagel@cern.ch>)  
**Module Type**: *GLOBAL*  
**Status**: Functional

### Description
The `Metronome` module can be used to slice data without strict event structure in arbitrarily long time frames, which serve as event definition for subsequent modules.
A new Event object is created and stored on the clipboard, with the begin and end of the frame calculated from the `event_length` parameter.

The module also provides the option to add trigger IDs to the generated event by specifying the number of triggers to be set per event via the `trigger`.
Trigger IDs are consecutive numbers, starting at zero.
With a setting of `triggers = 2`, the first event would receive trigger IDs 0 and 1, the subsequent event 2 and 3 and so forth.
A more detailed description is provided in the event building chapter of the user manual.

### Parameters
* `event_length`: Length of the event to be defined in physical units (not clock cycles of a specific device). Default value is `10us`.
* `skip_time`: Time to skip at the begin of the run before processing the first event. Defaults to `0us`.
* `triggers`: Number of triggers to generate and add to each event. All trigger timestamps are set to the center of the configured metronome time frame. Defaults to zero, i.e. no triggers are added.
* `skip_triggers`: Number of trigger IDs to skip at the begin of the run, before processing the first event. Defaults to 0. 

### Usage
```toml
[Metronome]
event_length = 500ns
```
