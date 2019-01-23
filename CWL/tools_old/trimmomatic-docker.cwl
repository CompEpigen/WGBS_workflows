dockerPull: dukegcb/trimmomatic
dockerFile: "#################################################################\n#\
  \ Dockerfile\n#\n# Software:         trimmomatic\n# Software Version: 0.32+dfsg-1\n\
  # Description:      DukeGCB trimmomatic image\n# Website:          http://www.usadellab.org/cms/?page=trimmomatic\n\
  # Provides:         trimmomatic\n# Base Image:       dukegcb/trimmomatic\n# Build\
  \ Cmd:        docker build --rm -t dukegcb/trimmomatic .\n# Pull Cmd:         docker\
  \ pull dukegcb/trimmomatic\n# Run Cmd:          docker run --rm -it dukegcb/trimmomatic\n\
  #################################################################\n\nFROM phusion/baseimage\n\
  MAINTAINER Dan Leehr <dan.leehr@duke.edu>\n\nRUN apt-get update && apt-get install\
  \ -y \\\n  openjdk-7-jre-headless \\\n  trimmomatic=\"0.32+dfsg-1\"\n\nCMD [\"/usr/bin/java\"\
  , \"-jar\", \"/usr/share/java/trimmomatic.jar\"]"
class: DockerRequirement

