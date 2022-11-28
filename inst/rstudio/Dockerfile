FROM jlaria/glasp_rstudio:latest
EXPOSE 8787
ENV PASSWORD=rstudio

COPY entrypoint.sh /usr/bin/
RUN chmod +x /usr/bin/entrypoint.sh
ENTRYPOINT ["entrypoint.sh"]
