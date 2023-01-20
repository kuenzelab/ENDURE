import streamlit as st
import streamlit_authenticator as stauth
# Import yaml
import yaml
# Import SafeLoader
from yaml import SafeLoader


#st.set_page_config(page_title="Settings", page_icon="🛠️")
with open('./config.yaml') as file:
    config = yaml.load(file, Loader=SafeLoader)
authenticator = stauth.Authenticate(
    config['credentials'],
    config['cookie']['name'],
    config['cookie']['key'],
    config['cookie']['expiry_days'],
    config['preauthorized']
    )
st.title("Settings")
# If user or username not in session state then set to None
if "name" not in st.session_state:
    st.session_state["name"] = None
if "username" not in st.session_state:
    st.session_state["username"] = None

# if the user is None then ask them to login
if st.session_state["name"] is None:
    # Please login
    st.write("Please go to the homepage to login")
else:
    # Create tabs for user settings
    tabs = [ "User Profile", "Password Reset"]
    tab1, tab2= st.tabs(tabs)



    with tab1:
        # Write the user's name to the session state
        st.write(f'Hello {st.session_state["name"]}!')
        st.write(f'Your username is {st.session_state["username"]}!')
        


    with tab2:
        st.write("Password Reset")
        if st.session_state["authentication_status"]:
            st.write(st.session_state["authentication_status"])
            try:
                #name, authentication_status, username = authenticator.login('Login', 'sidebar')
                if authenticator.reset_password(st.session_state["username"], 'Reset password'):
                    st.success('Password modified successfully')
            except Exception as e:
                st.error(e)

