#!/bin/bash
matlab -nodisplay -nosoftwarweopengl -nosplash -r \
  "NMSEvsSNR; EARvsNumberofMeasurements; "
