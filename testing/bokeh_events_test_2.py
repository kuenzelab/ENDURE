import streamlit as st
from bokeh.models import ColumnDataSource, CustomJS, TapTool, ResetTool,\
    LassoSelectTool
from bokeh.plotting import figure
import pandas as pd
import numpy as np


@st.cache
def data():
    info = pd.DataFrame(
        {
            "x": np.random.rand(500),
            "y": np.random.rand(500),
            "size": np.random.rand(500) * 10
        }
    )
    return info


def main():
    df = data()
    source = ColumnDataSource(df)
    st.subheader("Select Points From Map")
    callback = CustomJS(
        args=dict(source=source),
        code="""
        var select = source.selected.indices;
        console.log(select);
        """
    )

    plot = figure(tools=[ResetTool(), TapTool(), LassoSelectTool()])
    plot.circle(x="x", y="y", size="size", source=source, alpha=0.6)

    source.selected.js_on_change('indices', callback)
    st.bokeh_chart(plot)
    st.button('Hello World')


if __name__ == '__main__':
    main()
