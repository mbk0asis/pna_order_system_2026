# How to Deploy to Streamlit Community Cloud

The easiest way to deploy this app for free is using **Streamlit Community Cloud**.

## 1. Prepare your GitHub Repository
1.  Initialize Git in this folder (if not done):
    ```bash
    git init
    git add .
    git commit -m "Initial commit of PNA Designer"
    ```
2.  Create a **new repository** on [GitHub.com](https://github.com/new).
3.  Push your code:
    ```bash
    git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
    git branch -M main
    git push -u origin main
    ```

## 2. Deploy on Streamlit
1.  Go to [share.streamlit.io](https://share.streamlit.io/).
2.  Click **"New app"**.
3.  Select your GitHub repository (`YOUR_REPO_NAME`).
4.  Set **Main file path** to `app.py`.
5.  Click **"Deploy!"**.

## 3. Configuration (Optional)
If you want to hide development warnings or config options:
1.  In your Streamlit dashboard, go to "App Settings" > "Secrets".
2.  You can add any API keys if we add them later (currently not needed).

---

## Alternative: Deploy using Docker
If you prefer Docker, a `Dockerfile` can be created to containerize the application for deployment on AWS, Azure, or GCP.
