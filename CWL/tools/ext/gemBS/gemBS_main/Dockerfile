FROM heathsc/gembs
RUN rm start.sh
RUN echo "cd /opt; gemBS prepare -c example.conf -t example.csv \$@" > /home/start.sh
RUN echo "gemBS index" >> /home/start.sh
RUN echo "gemBS map" >> /home/start.sh
RUN echo "gemBS call" >> /home/start.sh
RUN echo "gemBS extract" >> /home/start.sh
RUN echo "gemBS map-report" >> /home/start.sh
RUN echo "gemBS call-report" >> /home/start.sh
ENTRYPOINT ["/bin/bash", "/home/start.sh"]
