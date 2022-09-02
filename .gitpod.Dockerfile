FROM gitpod/workspace-full

USER gitpod

# Install Julia
RUN pip install jill
RUN jill install --confirm

# Give control back to Gitpod Layer
USER root
