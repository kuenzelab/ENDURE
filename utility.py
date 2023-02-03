def load_text(folder_name: str, file_name: str) -> str:
    """
    Load text from file
    :param folder_name:
    :param file_name:
    :return:
    """
    with open(f'text/{folder_name}/{file_name}.txt', 'r') as file:
        data = file.read()
    return data.strip('\n')

import base64
import streamlit as st

@st.cache(allow_output_mutation=True)
def get_base64_of_bin_file(png_file):
    with open(png_file, "rb") as f:
        data = f.read()
    return base64.b64encode(data).decode()


def build_markup_for_logo(
    png_file,
    background_position="50% 10%",
    margin_top="10%",
    image_width="60%",
    image_height="200",
):
    binary_string = get_base64_of_bin_file(png_file)
    return """
            <style>
                [data-testid="stSidebarNav"] {
                    background-image: url("data:image/png;base64,%s");
                    background-repeat: no-repeat;
                    background-position: %s;
                    margin-top: %s;
                    background-size: %s %s;
                    padding-top: 60px;
                }
                [data-testid="stSidebarNav"]::before {
                content: "ENDURE";
                margin-left: 100px;
                margin-top: 20px;
                font-size: 30px;
                position: relative;
                top: 100px;
            }
            </style>
            """ % (
        binary_string,
        background_position,
        margin_top,
        image_width,
        image_height,
    )


def add_logo(png_file):
    logo_markup = build_markup_for_logo(png_file)
    st.markdown(
        logo_markup,
        unsafe_allow_html=True,
    )